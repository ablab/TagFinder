<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:variable name="width">500</xsl:variable>

    <xsl:template match="mzscan">
        <html>
            <title>Mzscan #<xsl:value-of select="scan-id"/></title>
            <body>
                <h3>Mzscan #<xsl:value-of select="scan-id"/></h3>


                <div>Precursor mass <xsl:value-of select="precursor-mass"/></div>


                <script>
                    var precursorMass = <xsl:value-of select="precursor-mass"/>;
                    var proton = 1.007276;
                    var width = 100;
                    function drawScale(ctx) {
                        ctx.strokeStyle = "red";
                        ctx.fillStyle = "red";
                        for (var i = -15; 35 >=i; i++) {
                            ctx.beginPath();
                            x = 500 + i * width;
                            ctx.moveTo(x, 500);
                            ctx.lineTo(x, 450);
                            ctx.stroke();
                            ctx.fillText(i, x +3, 470);
                        }
                    }
                    var mzpeaks = [<xsl:apply-templates select="mzpeak"/>];
                    var all = [];
                    var t = [];
                    for (var i = 0; mzpeaks.length > i; i++) {
                        t[i] = [];
                        for (var j = 1; 31 > j; j++) {
                            var m = mzpeaks[i][0] * j - proton * j;
                            if (m > 0 &amp;&amp; precursorMass + 100 > m) {
                                t[i][j] = {};
                                t[i][j].mass = m;
                                t[i][j].intencity = mzpeaks[i][1];
                                t[i][j].mzpeak = i;
                                t[i][j].charge = j;
                                t[i][j].score = 10;
                                all[all.length] = t[i][j];
                            }
                        }
                    }

                    function compareM(a,b) {
                        var am = a.mass;
                        var bm = b.mass;
                        if (bm > am) {
                            return -1;
                        }
                        if (am > bm) {
                            return 1;
                        }
                        return 0;
                    }
                    all.sort(compareM);
                </script>




                    <br/>
                    Scale <input type="range" id="scale" min="0" max="3000" value="100" style="width:500;" onchange="drawPeaks()"/><span id="scaleDisplay"></span>3000%
                    <a href="#" onclick="document.getElementById('scale').value=100;">Reset</a>
                    <br/>


                <canvas id="canvas" width="1000px" height="500"/>
<script>
    var ctx = document.getElementById("canvas").getContext('2d');

    var target;
    var emass;
    var modification;

    var pn = 1;
    var z = 1;



    function next() {
        var cur = 0;
        while (all.length > cur) {
            cur++;
            if (all[cur].mass > target) {
                if (all[cur].score > 1) {
                    break;
                }
            }
        }
        drawPeaks(all[cur].mzpeak, all[cur].charge);
    }

    function drawPeaks(pnParam, zParam) {
        var scale = document.getElementById('scale').value;
        if (pnParam != null) {
            pn = pnParam;
            z = zParam;
            target = t[pn][z].mass;
        }

        ctx.fillStyle = "white";
        ctx.clearRect(0, 0, 100000, 100000);
        ctx.strokeStyle = "black";
        ctx.fillStyle = "black";
        ctx.fillText("monomass : " + target, 800, 50);
        ctx.fillText("scale : " + scale + "%", 825, 60);
        drawScale(ctx);


        ctx.strokeStyle = "black";
        ctx.fillStyle = "black";

        var count = 0;
        for (var i = 0; mzpeaks.length > i; i++) {
            for (var j = 1; 31 > j; j++) {
                var p = t[i][j];
                if (p != null) {
                    if (p.score > 1) {
                        var m = p.mass;
                        var delta = m - target;
                        //if (10 > count) { alert(Math.abs(delta)); count++;}
                        if (15 > Math.abs(delta)) {
                            if (pn == i) {
                                ctx.strokeStyle = "green";
                                ctx.fillStyle = "green";
                            } else {
                                ctx.strokeStyle = "black";
                                ctx.fillStyle = "black";
                            }
                            ctx.beginPath();
                            x = 501 + (delta) * width;
                            ctx.moveTo(x, 500);
                            var y = 500 - p.intencity / 1000 * scale / 100;
                            ctx.lineTo(x, y);
                            ctx.stroke();
                            ctx.fillText(p.charge + " "  + delta.toFixed(2), x, y);
                        }
                    }
                }
            }
        }
    }

    function selectCurrentPeak() {
        for (var j = 1; 31 > j; j++) {
            if (t[pn][j] != null) {
                t[pn][j].score = j == z ? 100 : 0;
            }
        }
        for (var j = 1; 31 > j; j++) {
            if (t[pn][j] != null) {
                if (1 >= t[pn][j].score) {
                    table.rows[pn].cells[j - 1].innerHTML = "-";
                }
            }
        }
    
    }
</script>
                <br/>
                <a href="#" onclick="selectCurrentPeak(); return false;">select peak</a>
                <a href="#" onclick="next(); return false;">next</a>
                <br/>
                    <table id="table">
                    <!--<xsl:apply-templates select="mzpeak" mode="table"/>-->
                    </table>
                    <script>
                        var table = document.getElementById("table");
                        for (var i = 0; mzpeaks.length > i; i++) {
                            var row = table.insertRow(i);
                            for (var j = 1; 31 > j; j++) {
                                var cell = row.insertCell(j - 1);
                                var div = document.createElement('div')
                                div.innerHTML = "X";
                                if (t[i][j] != null) {
                                    div.innerHTML = '<a href="#" onclick="drawPeaks(' + i + ', ' + j + '); return false">' + t[i][j].mass + '</a>';
                                }
                                cell.appendChild(div);

                            }
                        }
                    </script>

            </body>
        </html>
    </xsl:template>


    <xsl:template match="mzpeak">
        [<xsl:value-of select="mass"/>, <xsl:value-of select="intencity"/>],
    </xsl:template>

    <xsl:template match="mzpeak" mode="table">
        <tr>
        <td><xsl:value-of select="mass"/></td>
        </tr>
    </xsl:template>


</xsl:stylesheet>