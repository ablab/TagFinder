function clear(ctx) {
    ctx.fillStyle = "white";
    ctx.clearRect(0, 0, 100000, 100000);
    ctx.strokeStyle = "black";
    ctx.fillStyle = "black";
}

var dist = [];
var startStat = [];
var finishStat = [];
var totalStat = [];
var cleavages = [];
var match = [];

var ppm = 10/1000/1000;

var deltas = [];

var modes = [];

function matchTen(peak, mass) {
     var error = mass * ppm;
     return peak.mass > mass - error && peak.mass < mass + error;
}

function matchRecalibrated(peak, mass, diff) {
     if (diff == null) {
        diff = 0;
     }
     var error = mass * ppm * 0.5;
     var m = peak.recalibrationMass - diff;
     return m > mass - error && m  < mass + error;
}

function initPrsm() {
    for (var i = 0; i <= sequence.length; i++) {
        dist[i] = [];
    }
    for (var i = 0; i <= sequence.length; i++) {
        cleavages[i] = i;
        startStat[i] = 0;
        finishStat[i] = 0;
        totalStat[i] = 0;

        var mass = 0;
        dist[i][i] = 0;
        for (var j = i; j < sequence.length; j++) {
            mass += getMass(sequence[j]);
            dist[i][j + 1] = mass;
            dist[j + 1][i] = mass;
            for (var p = 0; p < peaks.length; p++) {
                if (matchTen(peaks[p], mass)) {
                    deltas[deltas.length] = [peaks[p].mass - mass, peaks[p].mass, "\n"];
                }
            }
        }
    }

    scaleControl.value = Math.round((window.innerWidth - 100) * 100/dist[0][sequence.length]);

    var recalibration = 0;

    if (deltas.length >= 10) {
        deltas.sort(compareDeltas);
        limit = Math.floor(deltas.length / 10);
        var total = 0;
        for (var i = limit; i < deltas.length - limit; i++) {
            total += getErrorCoef(deltas[i]);
        }
        recalibration = total / (deltas.length - 2 * limit)
    }
    for (var p = 0; p < peaks.length; p++) {
        peaks[p].recalibrationMass = peaks[p].mass * (1 - recalibration);
    }

    deltas = [];
    var done = [];
    for (var p = 0; p < peaks.length; p++) {
        done[p] = 0;
    }
    for (var i = 0; i < sequence.length; i++) {
        match[i] = [];
        for (var j = i; j < sequence.length; j++) {
            match[i][j + 1] = 0;
            for (var p = 0; p < peaks.length; p++) {
                if (matchRecalibrated(peaks[p], dist[i][j + 1])) {
                    done[p]++;
                    deltas[deltas.length] = [peaks[p].recalibrationMass - dist[i][j + 1], peaks[p].recalibrationMass, "\n"];
                    if (match[i][j + 1] == 0) {
                        startStat[i]++;
                        finishStat[j + 1]++;
                        totalStat[i]++;
                        totalStat[j + 1]++;
                    }
                    match[i][j + 1]++;
                }
            }
        }
    }

    cleavages.sort(compareCleavages);
    var unmatched = [];
    for (var p = 0; p < peaks.length; p++) {
        if (done[p] == 0) {
            unmatched[unmatched.length] = peaks[p];
        }
    }


    do {
        bestScore = getBestScore(unmatched);
        if (bestScore.score < 2) {
            break;
        }
        modes[modes.length] = bestScore;
        unmatched = [];
        for (var u = 0; u < bestScore.unmatched.length; u++) {
            var peak = bestScore.unmatched[u];
            var described = 0;
            for (var j = 0; j <= sequence.length; j++) {
                if (matchRecalibrated(peak, dist[bestScore.cleavage][j], bestScore.diff)) {
                    described = 1;
                    break;
                }
            }
            if (described == 0) {
                unmatched[unmatched.length] = peak;
            }
        }
    } while (unmatched.length > 1);
}

function getBestScore(unmatched) {
    var bestScore = {diff:0, score:0};
    for (var i = 0; i <= sequence.length; i++) {
        var diff = [];
        for (var j = 0; j  <= sequence.length; j++) {
            for (var p = 0; p < unmatched.length; p++) {
                var peak = unmatched[p];
                diff[diff.length] = peak.recalibrationMass - dist[i][j];
            }
        }
        var shifts = getPossibleShifts(diff);
        for (var s = 0; s < shifts.length; s++) {
            var shift = shifts[s];
            var score = 0;
            var delta = 0;
            for (var u = 0; u < unmatched.length; u++) {
                var peak = unmatched[u];
                for (var j = 0; j <= sequence.length; j++) {
                    if (matchRecalibrated(peak, dist[i][j], shift)) {
                        score++;
                        delta += peak.recalibrationMass - shift - dist[i][j];
                        break;
                    }
                }
            }
            shift += delta / score;

            if (score > bestScore.score || (score == bestScore.score && Math.abs(shift) < Math.abs(bestScore.diff))) {
                bestScore.diff = shift;
                bestScore.score = score;
                bestScore.cleavage = i;
            }
        }
    }
    bestScore.unmatched = unmatched;
    return bestScore;
}


function getPossibleShifts(d) {
    d.sort(function(a,b){return a - b});
    var score = 0;
    var diff = 0;
    var temp = [];
    temp[0] = -1000000;
    var ans = [];

    for (var i = 1; i < d.length; i++) {
        if (d[i] - d[i - 1] < 1) {
            if (Math.abs(d[i]) < 100) {
                ans[ans.length] = (d[i] + d[i - 1])/2;
            }
        }
    }
    return ans;
}


function compareCleavages(c1, c2) {
    var ans = totalStat[c2] - totalStat[c1]
    if (ans == 0) {
        return c1 - c2;
    }
    return ans;
}

function compareDeltas(d1, d2) {
    return getErrorCoef(d1) - getErrorCoef(d2);
}

function getErrorCoef(d) {
    return d[0]/d[1];
}

var y;
var dx;
var k;

function repaintPrsm(scale) {
    clear(ctx);
    var mass = 0;
    if (scale <=0 || isNaN(scale)) {
        scale = 100;
    }
    k = scale /100;
    dx = dist[0][prefixLen];

    for (var i = 0; i <= sequence.length; i++) {
        ctx.beginPath();
        ctx.moveTo(k * (mass  - dx), 10);
        ctx.lineTo(k * (mass  - dx), 30);
        ctx.stroke();
        if (i < sequence.length) {
            var ch = sequence[i];
            var nextMass = mass + getMass(ch);
            ctx.fillText(ch, k * ((mass + nextMass)/2 - dx), 25);
            if (totalStat[i] > 0) {
                ctx.fillText(totalStat[i], k * (mass - dx), 45);
            }
            mass = nextMass;
        }
    }

    y = 47;
    ctx.strokeStyle = "green";
    //ctx.strokeStyle.lineWidth = 0.01;
    //alert(cleavages);
    for (var c = 0; totalStat[cleavages[c]] > 1; c++) {
        var i = cleavages[c];
        var yStart = y;
    //for (var i = 0; i < sequence.length; i++) {
        for (var j = 0; j < i; j++) {
            var count = 0;
            for (var p = 0; p < peaks.length; p++) {
                if (matchRecalibrated(peaks[p], dist[j][i])) {
                    count++;
                    drawLine(j, i, count);
                }
            }
        }

        for (var j = i + 1; j <= sequence.length; j++) {
            var count = 0;
            for (var p = 0; p < peaks.length; p++) {
                if (matchRecalibrated(peaks[p], dist[i][j])) {
                    count++;
                    drawLine(i, j, count);
                }
            }
        }
        ctx.save();
        ctx.strokeStyle = "red";
        ctx.beginPath();
        ctx.moveTo(k * (dist[0][i]  - dx), yStart + 2);
        ctx.lineTo(k * (dist[0][i]  - dx), y + 1);
        ctx.stroke();

        ctx.restore();
        y += 2;
    }

    for (var m = 0; m < modes.length; m++) {
        y += 5;
        var mod = modes[m];
        var i = mod.cleavage;
        var unmatched = mod.unmatched;
        ctx.fillText(mod.diff + " " + mod.score, k * (dist[0][i] - dx), y + 12);
        for (var j = 0; j <= sequence.length; j++) {
            var count = 0;
            for (var p = 0; p < unmatched.length; p++) {
                if (matchRecalibrated(unmatched[p], dist[i][j], mod.diff)) {
                    count++;
                    drawLine(i, j, count);
                }
            }
        }
    }
}

function drawLine(i, j, count) {
    if (j < i) {
        j += i;
        i = j - i;
        j -= i;
    }
    ctx.beginPath();
    if (count == 1) {
        y += 3;
        ctx.moveTo(k * (dist[0][i]  - dx), y);
        ctx.lineTo(k * (dist[0][j]  - dx), y);
        ctx.stroke();
    } else {
        ctx.arc(k * (dist[0][j]  - dx) - 5 * count, y - 1, 2, 0, Math.PI*2, false);
        ctx.closePath();
        ctx.fill();
    }
}

function getMass(ch) {
    switch (ch) {
        case 'G' : return 57.02146;
        case 'A' : return 71.03711;
        case 'S' : return 87.03203;
        case 'P' : return 97.05276;
        case 'V' : return 99.06841;
        case 'T' : return 101.04768;
        case 'C' : return 103.00919 + 57.021464; //CamC is a modification of cysteine
        case 'I' : return 113.08406;
        case 'L' : return 113.08406;
        case 'N' : return 114.04293;
        case 'D' : return 115.02694;
        case 'R' : return 156.10111;
        case 'Q' : return 128.05858;
        case 'K' : return 128.09496;
        case 'E' : return 129.04259;
        case 'M' : return 131.04049;
        case 'H' : return 137.05891;
        case 'F' : return 147.06841;
        case 'U' : return 150.953636;  //too rare
        case 'Y' : return 163.06333;
        case 'W' : return 186.07931;
        case 'O' : return 237.147727; //too rare
    }
    return 0;
}
