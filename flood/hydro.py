"""Depression filling, D8 flow routing, flow accumulation, carving and height above the drainage."""
from __future__ import annotations

import heapq
from dataclasses import dataclass

import numpy as np
from numba import njit

ROW_STEP = np.array([-1, -1, -1, 0, 0, 1, 1, 1])
COL_STEP = np.array([-1, 0, 1, -1, 1, -1, 0, 1])
STEP_LENGTH = np.array([2**0.5, 1.0, 2**0.5, 1.0, 1.0, 2**0.5, 1.0, 2**0.5])


@njit(cache=True)
def _priority_flood(z, closed_edge):
    """Fill depressions from the edge inwards (Barnes et al., 2014).

    Returns the filled surface, the cells in the order in which they were reached, and for each cell
    the cell from which it was reached. That order never goes down in filled elevation. The flood
    does not start from edge cells marked in `closed_edge`, so no flow leaves through them.
    """
    rows, cols = z.shape
    filled = z.copy()
    parent = np.full(rows * cols, -1, np.int64)
    closed = np.zeros((rows, cols), np.bool_)
    order = np.empty(rows * cols, np.int64)
    heap = [(0.0, 0, 0)]
    heap.pop()
    pushed = 0
    for r in range(rows):
        for c in range(cols):
            if np.isnan(z[r, c]):
                continue
            edge = r == 0 or c == 0 or r == rows - 1 or c == cols - 1
            if not edge:
                for k in range(8):
                    if np.isnan(z[r + ROW_STEP[k], c + COL_STEP[k]]):
                        edge = True
                        break
            if edge and not closed_edge[r, c]:
                closed[r, c] = True
                heapq.heappush(heap, (z[r, c], pushed, r * cols + c))
                pushed += 1
    count = 0
    while len(heap) > 0:
        level, _, index = heapq.heappop(heap)
        order[count] = index
        count += 1
        r, c = index // cols, index % cols
        for k in range(8):
            rr, cc = r + ROW_STEP[k], c + COL_STEP[k]
            if rr < 0 or cc < 0 or rr >= rows or cc >= cols or closed[rr, cc] or np.isnan(z[rr, cc]):
                continue
            closed[rr, cc] = True
            if filled[rr, cc] < level:
                filled[rr, cc] = level
            parent[rr * cols + cc] = index
            # the counter makes cells of equal level leave in the order they entered, which drains flats outwards
            heapq.heappush(heap, (filled[rr, cc], pushed, rr * cols + cc))
            pushed += 1
    return filled, order[:count], parent


@njit(cache=True)
def _receivers(filled, parent):
    """Steepest-descent neighbour of every cell; on flats, the cell from which the flood reached it."""
    rows, cols = filled.shape
    receiver = parent.copy()
    for r in range(rows):
        for c in range(cols):
            here = filled[r, c]
            if np.isnan(here):
                continue
            steepest = 0.0
            for k in range(8):
                rr, cc = r + ROW_STEP[k], c + COL_STEP[k]
                if rr < 0 or cc < 0 or rr >= rows or cc >= cols:
                    continue
                there = filled[rr, cc]
                if np.isnan(there):
                    continue
                gradient = (here - there) / STEP_LENGTH[k]
                if gradient > steepest:
                    steepest = gradient
                    receiver[r * cols + c] = rr * cols + cc
    return receiver


@njit(cache=True)
def _accumulate(order, receiver, weight):
    total = weight.copy()
    for position in range(len(order) - 1, -1, -1):
        index = order[position]
        if receiver[index] >= 0:
            total[receiver[index]] += total[index]
    return total


@njit(cache=True)
def _carve(order, receiver, z):
    """Lower every cell to the lowest elevation found upstream of it along the flow paths."""
    carved = z.copy()
    for position in range(len(order) - 1, -1, -1):
        index = order[position]
        target = receiver[index]
        if target >= 0 and carved[index] < carved[target]:
            carved[target] = carved[index]
    return carved


@njit(cache=True)
def _monotone_fit(order, receiver, z, stream):
    """Levels of the stream cells that never rise downstream and are closest to z by least squares.

    Isotonic regression on the stream tree by pooling. The cells are taken from the top of the streams
    down. A cell that lies above a block of cells upstream of it is pooled with the lowest such block,
    and the block takes the mean of its cells, until no block upstream is lower. For this order, in
    which every cell has one cell downstream, pooling the lowest violator gives the least-squares
    solution (Thompson, 1962). NaN outside the stream cells.
    """
    cells = order[stream[order]]                 # the stream cells, each after the cell it drains to
    count = len(cells)
    number = np.full(len(z), -1, np.int64)
    for k in range(count):
        number[cells[k]] = k
    total = np.empty(count)                      # sum of the elevations of the block headed by a cell
    size = np.ones(count)
    owner = np.full(count, -1, np.int64)         # the head of the block that took this block in
    first = np.full(count, -1, np.int64)         # blocks directly upstream of a block, as a linked list
    last = np.full(count, -1, np.int64)
    following = np.full(count, -1, np.int64)
    for head in range(count - 1, -1, -1):
        total[head] = z[cells[head]]
        while True:
            lowest, before_lowest, previous = -1, -1, -1
            block = first[head]
            while block >= 0:
                if lowest < 0 or total[block] * size[lowest] < total[lowest] * size[block]:
                    lowest, before_lowest = block, previous
                previous = block
                block = following[block]
            if lowest < 0 or total[lowest] * size[head] >= total[head] * size[lowest]:
                break
            if before_lowest < 0:
                first[head] = following[lowest]
            else:
                following[before_lowest] = following[lowest]
            if last[head] == lowest:
                last[head] = before_lowest
            if first[lowest] >= 0:
                if first[head] < 0:
                    first[head] = first[lowest]
                else:
                    following[last[head]] = first[lowest]
                last[head] = last[lowest]
            total[head] += total[lowest]
            size[head] += size[lowest]
            owner[lowest] = head
        target = receiver[cells[head]]
        if target >= 0 and stream[target]:
            below = number[target]
            following[head] = -1
            if first[below] < 0:
                first[below] = head
            else:
                following[last[below]] = head
            last[below] = head
    level = np.full(len(z), np.nan)
    value = np.empty(count)
    for k in range(count):
        value[k] = value[owner[k]] if owner[k] >= 0 else total[k] / size[k]
        level[cells[k]] = value[k]
    return level


@njit(cache=True)
def _value_downstream(order, receiver, value, stream):
    """For every cell, the value at the first stream cell on its flow path; NaN where there is none."""
    found = np.full(len(value), np.nan)
    for index in order:
        if stream[index]:
            found[index] = value[index]
        elif receiver[index] >= 0:
            found[index] = found[receiver[index]]
    return found


@dataclass
class FlowRouting:
    z: np.ndarray            # elevation as given, NaN outside the data
    filled: np.ndarray       # elevation with the depressions filled to their spill level
    carved: np.ndarray       # elevation with a channel cut through every barrier on a flow path
    order: np.ndarray        # flat cell indices, every cell after the cell it drains to
    receiver: np.ndarray     # flat index of the downstream neighbour, -1 where flow leaves the data
    cells: np.ndarray        # number of cells draining through each cell, the cell itself included

    def outlets(self):
        """Flat indices of the cells where flow leaves the data, largest contributing area first."""
        index = np.flatnonzero((self.receiver < 0) & np.isfinite(self.filled).ravel())
        return index[np.argsort(self.cells.ravel()[index])[::-1]]

    def flow_path(self, index):
        """Flat indices from a cell downstream to the cell where its flow leaves the data."""
        path = [int(index)]
        while self.receiver[path[-1]] >= 0:
            path.append(int(self.receiver[path[-1]]))
        return np.array(path)

    def step_length(self, index):
        """Distance from each of the given cells to its downstream neighbour, in cells (0 at an outlet)."""
        index = np.asarray(index)
        cols = self.filled.shape[1]
        target = self.receiver[index]
        length = np.hypot(index // cols - target // cols, index % cols - target % cols)
        return np.where(target >= 0, length, 0.0)


def route_flow(z, closed=None):
    """Route the flow of an elevation grid (NaN outside the data) with D8 on its depression-filled surface.

    The flow paths are those of the filled surface. The carved surface keeps the original elevations and
    lowers the cells of a barrier along each path instead (carving; Soille et al., 2003). Flow leaves the
    data at its edge, except through the cells marked in `closed`: mark the edge where a river enters.
    """
    z = np.ascontiguousarray(z, dtype=np.float64)
    closed = np.zeros(z.shape, np.bool_) if closed is None else np.ascontiguousarray(closed, dtype=np.bool_)
    filled, order, parent = _priority_flood(z, closed)
    if len(order) != np.isfinite(z).sum():
        raise ValueError("the closed edge leaves some cells without a way out of the data")
    receiver = _receivers(filled, parent)
    weight = np.isfinite(z).ravel().astype(np.float64)
    cells = _accumulate(order, receiver, weight).reshape(z.shape)
    cells[np.isnan(z)] = np.nan
    carved = _carve(order, receiver, z.ravel()).reshape(z.shape)
    return FlowRouting(z, filled, carved, order, receiver, cells)


def stream_level(routing, stream, level="fitted"):
    """Level of every stream cell by one of three readings of the elevations, NaN elsewhere.

    A stream must not rise downstream, and the elevations of its cells do. "filled" raises the low cells:
    it is the lowest level that never rises downstream and is nowhere below the elevations. "carved"
    lowers the high cells: it is the highest such level that is nowhere above them. "fitted" is the
    level that never rises downstream and is closest to the elevations of the stream cells by least
    squares; it lies between the other two.
    """
    stream = np.ascontiguousarray(stream, dtype=np.bool_) & np.isfinite(routing.z)
    if level == "filled":
        return np.where(stream, routing.filled, np.nan)
    if level == "carved":
        return np.where(stream, routing.carved, np.nan)
    if level == "fitted":
        return _monotone_fit(routing.order, routing.receiver, routing.z.ravel(), stream.ravel()).reshape(stream.shape)
    raise ValueError("level must be 'filled', 'fitted' or 'carved'")


def value_at_drainage(routing, stream, values):
    """For every cell, the value that `values` has at the stream cell the cell drains to.

    A stream cell takes its own value. NaN where the flow leaves the data without reaching a stream cell.
    """
    stream = np.ascontiguousarray(stream, dtype=np.bool_) & np.isfinite(routing.z)
    values = np.ascontiguousarray(values, dtype=np.float64).ravel()
    found = _value_downstream(routing.order, routing.receiver, values, stream.ravel())
    return found.reshape(routing.z.shape)


def height_above_drainage(routing, stream, level="fitted"):
    """Height of the ground above the level of the stream cell it drains to (Renno et al., 2008).

    The ground is the elevation as given; the level of the stream cells is that of `stream_level`.
    Stream cells are at height zero. The height is negative where the ground lies below the level given
    to its stream, and NaN where the flow leaves the data without reaching a stream cell.
    """
    stream = np.ascontiguousarray(stream, dtype=np.bool_) & np.isfinite(routing.z)
    height = routing.z - value_at_drainage(routing, stream, stream_level(routing, stream, level))
    height[stream] = 0.0
    return height
