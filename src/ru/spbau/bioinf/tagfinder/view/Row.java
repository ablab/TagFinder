package ru.spbau.bioinf.tagfinder.view;

import java.util.HashMap;

public class Row {
    private int row;
    private HashMap<Integer, Cell> cells = new HashMap<Integer, Cell>();

    public Row(int row) {
        this.row = row;
    }

    public HashMap<Integer, Cell> getCells() {
        return cells;
    }

    public Content getContent(int col) {
        if (cells.containsKey(col)) {
            return cells.get(col).getContent();
        }
        return null;
    }

    public void addCell(int col, Content content) {
        cells.put(col, new Cell(row, col, content));
    }
}
