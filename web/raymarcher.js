/*
 * WebGL Implicit Algebraic Surface Raymarcher
 * (c) 2026 Aneesh Krishna
 *
 * Implements Taubin first-order implicit distance estimation:
 * d(p) approx |f(p)| / ||grad f(p)||
 * for Whitney umbrella, Cayley cubic, Clebsch cubic, Kummer quartic,
 * Barth sextic, Enneper surface, and Calabi-Yau quintic cross sections.
 */

(function () {
  const canvas = document.getElementById("glCanvas");
  const gl = canvas.getContext("webgl");

  if (!gl) {
    alert("WebGL not supported by your browser.");
    return;
  }

  // Vertex Shader
  const vsSource = `
    attribute vec2 aPosition;
    void main() {
      gl_Position = vec4(aPosition, 0.0, 1.0);
    }
  `;

  // Fragment Shader
  const fsSource = `
    precision highp float;
    uniform vec2 uResolution;
    uniform vec2 uMouse;
    uniform float uTime;
    uniform float uParamC;
    uniform float uCameraDist;
    uniform int uSurfaceType;
    uniform int uColorMode;

    // Polynomial evaluations for implicit algebraic varieties V(f) = 0
    float evalSurface(vec3 p) {
      float x = p.x;
      float y = p.y;
      float z = p.z;
      float c = uParamC;

      if (uSurfaceType == 0) {
        // Whitney Umbrella: x^2 - y^2 * z = c
        return (x * x - y * y * z) - c;
      } else if (uSurfaceType == 1) {
        // Cayley Cubic: 4*(x^2 + y^2 + z^2) + 16*x*y*z - 1 = c
        return (4.0 * (x * x + y * y + z * z) + 16.0 * x * y * z - 1.0) - c;
      } else if (uSurfaceType == 2) {
        // Clebsch Diagonal Cubic: 81*(x^3 + y^3 + z^3) - 189*(x^2*y + ...) + 54*x*y*z + ...
        float s1 = x + y + z;
        float s2 = x*y + y*z + z*x;
        float s3 = x*y*z;
        return (81.0 * (s1*s1*s1 - 3.0*s1*s2 + 3.0*s3) - 189.0 * (s1*s2 - 3.0*s3) + 54.0 * s3 + 126.0 * s2 - 9.0 * (s1*s1 - 2.0*s2) - 9.0 * s1 + 1.0) - c;
      } else if (uSurfaceType == 3) {
        // Kummer Quartic: (x^2 + y^2 + z^2 - w^2)^2 - (3*w^2 - 1)/(3 - w^2) * ...
        float r2 = x*x + y*y + z*z;
        return (r2*r2 - 2.0*r2 + 0.3 + 4.0*x*y*z) - c;
      } else if (uSurfaceType == 4) {
        // Barth Sextic (simplified core polynomial): 4*(phi^2 * x^2 - y^2)*(phi^2 * y^2 - z^2)*(phi^2 * z^2 - x^2) - (1 + 2*phi)*(x^2 + y^2 + z^2 - 1)^2
        float phi = 1.6180339887;
        float f1 = (phi*phi*x*x - y*y);
        float f2 = (phi*phi*y*y - z*z);
        float f3 = (phi*phi*z*z - x*x);
        float r2 = x*x + y*y + z*z;
        return (4.0 * f1 * f2 * f3 - (1.0 + 2.0*phi) * (r2 - 1.0) * (r2 - 1.0)) - c;
      } else if (uSurfaceType == 5) {
        // Enneper Minimal Surface algebraic form
        float r2 = x*x + y*y + z*z;
        return (64.0 * z*z*z + 48.0 * z * (x*x + y*y) - (x*x + y*y)*(x*x + y*y)) - c;
      } else if (uSurfaceType == 6) {
        // Calabi-Yau Fermat Quintic slice: Re(z_1^5 + z_2^5) = c
        float x5 = x*x*x*x*x - 10.0*x*x*x*y*y + 5.0*x*y*y*y*y;
        float z5 = z*z*z*z*z;
        return (x5 + z5 - 0.2) - c;
      } else {
        // Smooth Torus (R = 1.2, r = 0.45)
        float R = 1.2;
        float r = 0.45;
        float d1 = sqrt(x*x + y*y) - R;
        return (d1*d1 + z*z - r*r) - c;
      }
    }

    // Finite difference gradient for Taubin distance estimator:
    // d(p) approx |f(p)| / ||grad f(p)||
    vec3 calcNormal(vec3 p) {
      float eps = 0.002;
      float f0 = evalSurface(p);
      float dx = evalSurface(vec3(p.x + eps, p.y, p.z)) - f0;
      float dy = evalSurface(vec3(p.x, p.y + eps, p.z)) - f0;
      float dz = evalSurface(vec3(p.x, p.y, p.z + eps)) - f0;
      return normalize(vec3(dx, dy, dz));
    }

    float taubinDistance(vec3 p) {
      float fVal = evalSurface(p);
      float eps = 0.003;
      float dx = (evalSurface(vec3(p.x + eps, p.y, p.z)) - evalSurface(vec3(p.x - eps, p.y, p.z))) / (2.0 * eps);
      float dy = (evalSurface(vec3(p.x, p.y + eps, p.z)) - evalSurface(vec3(p.x, p.y - eps, p.z))) / (2.0 * eps);
      float dz = (evalSurface(vec3(p.x, p.y, p.z + eps)) - evalSurface(vec3(p.x, p.y, p.z - eps))) / (2.0 * eps);
      float gradLen = sqrt(dx*dx + dy*dy + dz*dz + 1e-6);
      return abs(fVal) / gradLen;
    }

    void main() {
      vec2 uv = (gl_FragCoord.xy - 0.5 * uResolution.xy) / uResolution.y;

      // Orbit camera rotation
      float yaw = uMouse.x * 3.14159 * 2.0;
      float pitch = clamp(uMouse.y * 3.14159 - 1.57, -1.5, 1.5);

      vec3 ro = vec3(
        uCameraDist * cos(pitch) * sin(yaw),
        uCameraDist * sin(pitch),
        uCameraDist * cos(pitch) * cos(yaw)
      );
      vec3 target = vec3(0.0, 0.0, 0.0);
      vec3 ww = normalize(target - ro);
      vec3 uu = normalize(cross(ww, vec3(0.0, 1.0, 0.0)));
      vec3 vv = cross(uu, ww);
      vec3 rd = normalize(uv.x * uu + uv.y * vv + 1.6 * ww);

      // Sphere tracing loop
      float t = 0.1;
      float maxDist = 8.0;
      float hitDist = 0.003;
      bool hit = false;
      vec3 hitPos = ro;

      for (int i = 0; i < 96; ++i) {
        vec3 p = ro + t * rd;
        float d = taubinDistance(p);
        if (d < hitDist) {
          hit = true;
          hitPos = p;
          break;
        }
        t += max(0.002, 0.45 * d);
        if (t > maxDist) break;
      }

      if (hit) {
        vec3 norm = calcNormal(hitPos);
        vec3 lightDir = normalize(vec3(0.6, 0.9, 0.5));
        float diff = max(0.0, dot(norm, lightDir));
        float amb = 0.2;

        vec3 color = vec3(0.0);
        if (uColorMode == 0) {
          // Normal gradient
          color = norm * 0.5 + 0.5;
        } else if (uColorMode == 1) {
          // Studio diffuse
          vec3 baseColor = vec3(0.35, 0.65, 0.95);
          color = baseColor * (amb + 0.8 * diff) + vec3(0.3) * pow(max(0.0, dot(reflect(-lightDir, norm), -rd)), 16.0);
        } else if (uColorMode == 2) {
          // Curvature / gradient
          float cur = length(norm.xy);
          color = vec3(cur, 1.0 - cur, 0.8);
        } else {
          // Depth shading
          float depthVal = 1.0 - (t / maxDist);
          color = vec3(depthVal * 0.4, depthVal * 0.7, depthVal);
        }

        gl_FragColor = vec4(color, 1.0);
      } else {
        // Gradient background
        vec3 topColor = vec3(0.04, 0.05, 0.08);
        vec3 botColor = vec3(0.01, 0.02, 0.03);
        vec3 bg = mix(botColor, topColor, uv.y + 0.5);
        gl_FragColor = vec4(bg, 1.0);
      }
    }
  `;

  function createShader(gl, type, source) {
    const shader = gl.createShader(type);
    gl.shaderSource(shader, source);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
      console.error(gl.getShaderInfoLog(shader));
      gl.deleteShader(shader);
      return null;
    }
    return shader;
  }

  const vs = createShader(gl, gl.VERTEX_SHADER, vsSource);
  const fs = createShader(gl, gl.FRAGMENT_SHADER, fsSource);
  const program = gl.createProgram();
  gl.attachShader(program, vs);
  gl.attachShader(program, fs);
  gl.linkProgram(program);

  if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
    console.error(gl.getProgramInfoLog(program));
    return;
  }

  gl.useProgram(program);

  // Screen quad
  const quadBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, quadBuf);
  gl.bufferData(
    gl.ARRAY_BUFFER,
    new Float32Array([-1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1]),
    gl.STATIC_DRAW
  );

  const aPosLoc = gl.getAttribLocation(program, "aPosition");
  gl.enableVertexAttribArray(aPosLoc);
  gl.vertexAttribPointer(aPosLoc, 2, gl.FLOAT, false, 0, 0);

  // Uniform locations
  const uResLoc = gl.getUniformLocation(program, "uResolution");
  const uMouseLoc = gl.getUniformLocation(program, "uMouse");
  const uTimeLoc = gl.getUniformLocation(program, "uTime");
  const uParamCLoc = gl.getUniformLocation(program, "uParamC");
  const uCamDistLoc = gl.getUniformLocation(program, "uCameraDist");
  const uSurfTypeLoc = gl.getUniformLocation(program, "uSurfaceType");
  const uColorModeLoc = gl.getUniformLocation(program, "uColorMode");

  // Interaction State
  let mouse = [0.25, 0.35];
  let isDragging = false;
  let lastX = 0, lastY = 0;

  canvas.addEventListener("mousedown", (e) => {
    isDragging = true;
    lastX = e.clientX;
    lastY = e.clientY;
  });

  window.addEventListener("mouseup", () => {
    isDragging = false;
  });

  window.addEventListener("mousemove", (e) => {
    if (!isDragging) return;
    const dx = e.clientX - lastX;
    const dy = e.clientY - lastY;
    lastX = e.clientX;
    lastY = e.clientY;
    mouse[0] -= dx * 0.005;
    mouse[1] += dy * 0.005;
  });

  // UI inputs
  const surfaceSelect = document.getElementById("surfaceSelect");
  const paramC = document.getElementById("paramC");
  const valC = document.getElementById("valC");
  const paramZoom = document.getElementById("paramZoom");
  const valZoom = document.getElementById("valZoom");
  const colorModeSelect = document.getElementById("colorModeSelect");

  const statSurface = document.getElementById("statSurface");
  const statDegree = document.getElementById("statDegree");
  const statFPS = document.getElementById("statFPS");

  const surfaceNames = [
    ["Whitney Umbrella", "3"],
    ["Cayley Cubic", "3"],
    ["Clebsch Diagonal Cubic", "3"],
    ["Kummer Quartic", "4"],
    ["Barth Sextic", "6"],
    ["Enneper Minimal Surface", "9"],
    ["Calabi-Yau Quintic Slice", "5"],
    ["Smooth Torus", "4"]
  ];

  paramC.addEventListener("input", (e) => {
    valC.textContent = parseFloat(e.target.value).toFixed(2);
  });

  paramZoom.addEventListener("input", (e) => {
    valZoom.textContent = parseFloat(e.target.value).toFixed(1);
  });

  surfaceSelect.addEventListener("change", (e) => {
    const idx = parseInt(e.target.value, 10);
    statSurface.textContent = surfaceNames[idx][0];
    statDegree.textContent = surfaceNames[idx][1];
  });

  // Render Loop
  let lastTime = performance.now();
  let frameCount = 0;
  let fpsTimer = performance.now();

  function resize() {
    const w = canvas.clientWidth;
    const h = canvas.clientHeight;
    if (canvas.width !== w || canvas.height !== h) {
      canvas.width = w;
      canvas.height = h;
      gl.viewport(0, 0, w, h);
    }
  }

  function render(time) {
    resize();

    frameCount++;
    if (time - fpsTimer >= 1000) {
      statFPS.textContent = frameCount.toString();
      frameCount = 0;
      fpsTimer = time;
    }

    gl.uniform2f(uResLoc, canvas.width, canvas.height);
    gl.uniform2f(uMouseLoc, mouse[0], mouse[1]);
    gl.uniform1f(uTimeLoc, time * 0.001);
    gl.uniform1f(uParamCLoc, parseFloat(paramC.value));
    gl.uniform1f(uCamDistLoc, parseFloat(paramZoom.value));
    gl.uniform1i(uSurfTypeLoc, parseInt(surfaceSelect.value, 10));
    gl.uniform1i(uColorModeLoc, parseInt(colorModeSelect.value, 10));

    gl.drawArrays(gl.TRIANGLES, 0, 6);
    requestAnimationFrame(render);
  }

  requestAnimationFrame(render);
})();
