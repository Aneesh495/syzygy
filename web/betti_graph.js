/*
 * Graded Free Resolution Chain Complex DAG Visualizer
 * (c) 2026 Aneesh Krishna
 *
 * Renders the sequence of free modules and differential maps:
 * 0 <- M <- F_0 <- F_1 <- F_2 <- ... <- F_p <- 0
 * with generator degree shifts, Betti numbers, and syzygy ranks.
 */

(function () {
  const canvas = document.getElementById("graphCanvas");
  if (!canvas) return;
  const ctx = canvas.getContext("2d");

  // Free resolution data for preset surfaces
  const resolutions = [
    {
      // Whitney umbrella: x^2 - y^2*z = 0 (hypersurface in C^3)
      ranks: [1, 1],
      labels: ["R", "R(-3)"],
      betti: "       total:    1    1\n       -----------------\n           0:    1    .\n           1:    .\n           2:    .    1"
    },
    {
      // Twisted cubic curve in P^3: I = <y^2 - x*z, x*w - y*z, z^2 - y*w>
      ranks: [1, 3, 2],
      labels: ["R", "R^3(-2)", "R^2(-3)"],
      betti: "       total:    1    3    2\n       ---------------------\n           0:    1    .    .\n           1:    .    3    2"
    },
    {
      // Complete intersection of two quadrics: I = <q_1, q_2>
      ranks: [1, 2, 1],
      labels: ["R", "R^2(-2)", "R(-4)"],
      betti: "       total:    1    2    1\n       ---------------------\n           0:    1    .    .\n           1:    .    2    .\n           2:    .    .    1"
    }
  ];

  let currentRes = resolutions[1]; // default twisted cubic

  function resize() {
    const w = canvas.clientWidth;
    const h = canvas.clientHeight;
    if (canvas.width !== w || canvas.height !== h) {
      canvas.width = w;
      canvas.height = h;
    }
  }

  function drawGraph() {
    resize();
    const w = canvas.width;
    const h = canvas.height;
    ctx.clearRect(0, 0, w, h);

    const ranks = currentRes.ranks;
    const n = ranks.length;
    const padding = 40;
    const spacing = (w - 2 * padding) / (n + 1);

    // Draw connecting arrows
    for (let i = 0; i < n; ++i) {
      const x1 = padding + (i + 1) * spacing;
      const x2 = padding + i * spacing;
      const y = h / 2;

      ctx.strokeStyle = "#58a6ff";
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.moveTo(x1 - 24, y);
      ctx.lineTo(x2 + 24, y);
      ctx.stroke();

      // Arrowhead pointing left towards 0
      ctx.fillStyle = "#58a6ff";
      ctx.beginPath();
      ctx.moveTo(x2 + 24, y);
      ctx.lineTo(x2 + 32, y - 5);
      ctx.lineTo(x2 + 32, y + 5);
      ctx.closePath();
      ctx.fill();

      // Label differential d_i
      if (i > 0) {
        ctx.fillStyle = "#8b949e";
        ctx.font = "10px monospace";
        ctx.textAlign = "center";
        ctx.fillText("d_" + i, (x1 + x2) / 2, y - 10);
      }
    }

    // Draw module nodes
    for (let i = 0; i < n; ++i) {
      const x = padding + (i + 1) * spacing;
      const y = h / 2;

      // Glow circle
      ctx.shadowColor = "#1f6feb";
      ctx.shadowBlur = 10;
      ctx.fillStyle = "#12161f";
      ctx.beginPath();
      ctx.arc(x, y, 20, 0, Math.PI * 2);
      ctx.fill();

      // Border
      ctx.shadowBlur = 0;
      ctx.strokeStyle = i === 0 ? "#3fb950" : "#bc8cff";
      ctx.lineWidth = 2;
      ctx.stroke();

      // Text label
      ctx.fillStyle = "#ffffff";
      ctx.font = "bold 11px monospace";
      ctx.textAlign = "center";
      ctx.textBaseline = "middle";
      ctx.fillText("F_" + i, x, y);

      // Degree label below node
      ctx.fillStyle = "#8b949e";
      ctx.font = "9px monospace";
      ctx.fillText(currentRes.labels[i], x, y + 30);
    }

    // Target quotient module M = R/I at the left
    const x0 = padding;
    const y0 = h / 2;
    ctx.strokeStyle = "#d29922";
    ctx.lineWidth = 2;
    ctx.strokeRect(x0 - 18, y0 - 14, 36, 28);

    ctx.fillStyle = "#d29922";
    ctx.font = "bold 10px monospace";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText("R/I", x0, y0);
  }

  // Update when surface changes
  const surfaceSelect = document.getElementById("surfaceSelect");
  const bettiTableOutput = document.getElementById("bettiTableOutput");

  if (surfaceSelect) {
    surfaceSelect.addEventListener("change", (e) => {
      const idx = parseInt(e.target.value, 10);
      if (idx === 0) {
        currentRes = resolutions[0];
      } else if (idx === 1 || idx === 2) {
        currentRes = resolutions[1];
      } else {
        currentRes = resolutions[2];
      }
      if (bettiTableOutput) {
        bettiTableOutput.textContent = currentRes.betti;
      }
      drawGraph();
    });
  }

  drawGraph();
  window.addEventListener("resize", drawGraph);
})();
