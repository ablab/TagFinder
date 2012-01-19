package ru.spbau.bioinf.palign;

import ru.spbau.bioinf.tagfinder.Consts;
import ru.spbau.bioinf.tagfinder.Peak;

public enum ModificationType {
    WATER_LOSS(-Consts.WATER), ONE_LOSS(-1), NONE(0), ONE_GAIN(+1);

    double delta;

    private ModificationType(double delta) {
        this.delta = delta;
    }
    
    public double getError(Peak peak, double peptideMass) {
        return peak.getMass() - delta - peptideMass;
    }
}
