function newCell = deleteCellContents(oldCell, indices)

oldCell(indices) = [];
newCell = oldCell;