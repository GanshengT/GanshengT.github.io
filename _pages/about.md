---
permalink: /
title: "I study how bodily rhythms structure the neural dynamics of human cognition."
excerpt: "About me"
author_profile: true
redirect_from: 
  - /about/
  - /about.html
---

<style>
  .about-rhythm {
    margin-bottom: 2.5rem;
  }

  .about-fold {
    margin: 1.25rem 0;
    border: 1px solid rgba(15, 23, 42, 0.08);
    border-radius: 18px;
    background: linear-gradient(180deg, #fffefb 0%, #f5f7fb 100%);
    box-shadow: 0 12px 28px rgba(15, 23, 42, 0.05);
    overflow: hidden;
  }

  .about-fold summary {
    padding: 1rem 1.2rem;
    cursor: pointer;
    font-weight: 700;
    list-style: none;
  }

  .about-fold summary::-webkit-details-marker {
    display: none;
  }

  .about-fold summary::after {
    content: "+";
    float: right;
    color: #1f6a73;
    font-size: 1.15rem;
    line-height: 1;
  }

  .about-fold[open] summary::after {
    content: "−";
  }

  .about-fold__content {
    padding: 0 1.2rem 1.2rem;
  }

  .about-rhythm__visual {
    position: relative;
    min-height: 260px;
    border-radius: 24px;
    background:
      radial-gradient(circle at top left, rgba(243, 181, 107, 0.18), transparent 30%),
      radial-gradient(circle at bottom right, rgba(70, 159, 178, 0.22), transparent 34%),
      linear-gradient(135deg, #0f172a 0%, #123449 52%, #1f6a73 100%);
    overflow: hidden;
    box-shadow: 0 24px 54px rgba(15, 23, 42, 0.14);
  }

  .about-rhythm__canvas {
    position: absolute;
    inset: 0;
    width: 100%;
    height: 100%;
  }

  .about-rhythm__overlay {
    position: relative;
    z-index: 1;
    min-height: 260px;
  }

  @media (max-width: 48em) {
    .about-rhythm__visual,
    .about-rhythm__overlay {
      min-height: 220px;
    }
  }
</style>

<div class="about-rhythm" markdown="0">
  <section class="about-rhythm__visual" aria-label="Animated oscillation signal">
    <canvas class="about-rhythm__canvas" id="about-rhythm-canvas"></canvas>
    <div class="about-rhythm__overlay"></div>
  </section>
</div>

<script>
  document.addEventListener("DOMContentLoaded", function() {
    const canvas = document.getElementById("about-rhythm-canvas");
    if (!canvas) return;

    const context = canvas.getContext("2d");
    const visual = canvas.parentElement;
    const state = {
      tracks: {},
      fixationIndex: 0,
      nextIndex: 1,
      phase: "fixation",
      phaseStartedAt: 0,
      fixationDuration: 1600,
      saccadeDuration: 340,
      trails: [],
      events: [],
      sequence: [0, 3, 1, 4, 2, 5],
      sampleSpacingPx: 2,
      sampleStepMs: 22,
      signalAccumulator: 0,
      lastSignalTime: 0,
      heartPhase: 0,
      neuralPhaseSlow: 0,
      neuralPhaseFast: 0,
      ecgBuffer: [],
      neuralBuffer: []
    };

    function resize() {
      const ratio = window.devicePixelRatio || 1;
      canvas.width = visual.clientWidth * ratio;
      canvas.height = visual.clientHeight * ratio;
      canvas.style.width = visual.clientWidth + "px";
      canvas.style.height = visual.clientHeight + "px";
      context.setTransform(ratio, 0, 0, ratio, 0, 0);
      buildLayout();
    }

    function buildLayout() {
      const width = visual.clientWidth;
      const height = visual.clientHeight;
      const left = 78;
      const right = width - 28;
      const gazeY = height * 0.24;
      const ecgY = height * 0.52;
      const neuralY = height * 0.78;
      const span = Math.max(right - left, 180);
      const step = span / 5;

      state.tracks = {
        width: width,
        height: height,
        left: left,
        right: right,
        gazeY: gazeY,
        ecgY: ecgY,
        neuralY: neuralY,
        sampleCount: Math.max(90, Math.floor((right - left) / state.sampleSpacingPx)),
        gazeTargets: [
          { x: left + step * 0.2, y: gazeY - 10 },
          { x: left + step * 1.1, y: gazeY + 8 },
          { x: left + step * 2.0, y: gazeY - 14 },
          { x: left + step * 2.9, y: gazeY + 6 },
          { x: left + step * 4.0, y: gazeY - 9 },
          { x: left + step * 4.8, y: gazeY + 10 }
        ]
      };

      resetSignalBuffers();
    }

    function lerp(start, end, progress) {
      return start + (end - start) * progress;
    }

    function currentTarget(index) {
      const safeIndex = state.sequence[index % state.sequence.length];
      return state.tracks.gazeTargets[safeIndex];
    }

    function updateGaze(time) {
      if (!state.phaseStartedAt) {
        state.phaseStartedAt = time;
      }

      const elapsed = time - state.phaseStartedAt;
      const from = currentTarget(state.fixationIndex);
      const to = currentTarget(state.nextIndex);

      if (state.phase === "fixation") {
        if (elapsed >= state.fixationDuration) {
          state.phase = "saccade";
          state.phaseStartedAt = time;
          state.events.push({ ageMs: 0 });
        }
        return from;
      }

      const progress = Math.min(elapsed / state.saccadeDuration, 1);
      const eased = 1 - Math.pow(1 - progress, 3);
      const position = {
        x: lerp(from.x, to.x, eased),
        y: lerp(from.y, to.y, eased)
      };

      if (progress >= 1) {
        state.trails.push({ from: from, to: to, bornAt: time });
        state.fixationIndex = state.nextIndex;
        state.nextIndex += 1;
        state.phase = "fixation";
        state.phaseStartedAt = time;
      }

      return position;
    }

    function drawTrackLabel(text, y) {
      context.fillStyle = "rgba(248, 250, 252, 0.62)";
      context.font = "12px sans-serif";
      context.textBaseline = "middle";
      context.fillText(text, 18, y);
    }

    function drawGuide(y) {
      context.strokeStyle = "rgba(255, 255, 255, 0.08)";
      context.lineWidth = 1;
      context.beginPath();
      context.moveTo(state.tracks.left, y);
      context.lineTo(state.tracks.right, y);
      context.stroke();
    }

    function drawGazeTrack(time, gazePosition) {
      drawTrackLabel("Gaze", state.tracks.gazeY);
      drawGuide(state.tracks.gazeY);

      state.tracks.gazeTargets.forEach(function(target) {
        context.beginPath();
        context.fillStyle = "rgba(244, 246, 248, 0.16)";
        context.arc(target.x, target.y, 4, 0, Math.PI * 2);
        context.fill();
      });

      state.trails = state.trails.filter(function(trail) {
        return time - trail.bornAt < 1200;
      });

      state.trails.forEach(function(trail) {
        const alpha = 1 - (time - trail.bornAt) / 1200;
        context.strokeStyle = "rgba(243, 181, 107, " + (alpha * 0.5).toFixed(3) + ")";
        context.lineWidth = 1.4;
        context.beginPath();
        context.moveTo(trail.from.x, trail.from.y);
        context.lineTo(trail.to.x, trail.to.y);
        context.stroke();
      });

      if (state.phase === "saccade") {
        const from = currentTarget(state.fixationIndex);
        context.strokeStyle = "rgba(243, 181, 107, 0.55)";
        context.lineWidth = 1.25;
        context.beginPath();
        context.moveTo(from.x, from.y);
        context.lineTo(gazePosition.x, gazePosition.y);
        context.stroke();
      }

      context.beginPath();
      context.fillStyle = "rgba(255, 248, 233, 0.96)";
      context.shadowColor = "rgba(243, 181, 107, 0.85)";
      context.shadowBlur = 14;
      context.arc(gazePosition.x, gazePosition.y, 5.5, 0, Math.PI * 2);
      context.fill();
      context.shadowBlur = 0;
    }

    function bump(value, center, width, amplitude) {
      const normalized = (value - center) / width;
      return amplitude * Math.exp(-0.5 * normalized * normalized);
    }

    function ecgSample(phase) {
      let value = 0;

      value += bump(phase, 0.18, 0.028, 3.6);   // P wave
      value += bump(phase, 0.39, 0.010, -4.5);  // Q wave
      value += bump(phase, 0.42, 0.006, 22.5);  // R wave
      value += bump(phase, 0.45, 0.012, -8.8);  // S wave
      value += bump(phase, 0.70, 0.060, 6.8);   // T wave

      return value;
    }

    function heartbeatCoupling(phase) {
      return bump(phase, 0.42, 0.028, 1.1) + bump(phase, 0.70, 0.09, 0.35);
    }

    function gazeEvokedResponse(ageMs) {
      if (ageMs < 0 || ageMs > 900) {
        return 0;
      }

      const earlyNegative = bump(ageMs, 70, 24, -4.5);
      const positivePeak = bump(ageMs, 145, 36, 8.5);
      const burstEnvelope = Math.exp(-Math.max(ageMs - 110, 0) / 180);
      const burst = ageMs > 80 ? Math.sin(ageMs * 0.105) * 4.6 * burstEnvelope : 0;

      return earlyNegative + positivePeak + burst;
    }

    function stepSignals(steps) {
      for (let i = 0; i < steps; i += 1) {
        const ecg = ecgSample(state.heartPhase);
        const coupling = heartbeatCoupling(state.heartPhase);
        let evoked = 0;

        state.events.forEach(function(event) {
          evoked += gazeEvokedResponse(event.ageMs);
          event.ageMs += state.sampleStepMs;
        });

        state.events = state.events.filter(function(event) {
          return event.ageMs <= 900;
        });

        const oscillation =
          Math.sin(state.neuralPhaseSlow) * 3.8 +
          Math.sin(state.neuralPhaseFast) * 2.1;
        const neural = oscillation * (1 + coupling) + ecg * 0.11 + evoked;

        state.ecgBuffer.unshift(ecg);
        state.neuralBuffer.unshift(neural);

        if (state.ecgBuffer.length > state.tracks.sampleCount) {
          state.ecgBuffer.pop();
        }

        if (state.neuralBuffer.length > state.tracks.sampleCount) {
          state.neuralBuffer.pop();
        }

        state.heartPhase = (state.heartPhase + state.sampleStepMs / 1080) % 1;
        state.neuralPhaseSlow += 0.16;
        state.neuralPhaseFast += 0.45;
      }
    }

    function resetSignalBuffers() {
      state.signalAccumulator = 0;
      state.lastSignalTime = 0;
      state.heartPhase = 0;
      state.neuralPhaseSlow = 0;
      state.neuralPhaseFast = 0;
      state.events = [];
      state.ecgBuffer = [];
      state.neuralBuffer = [];

      stepSignals(state.tracks.sampleCount);
    }

    function advanceSignals(time) {
      if (!state.lastSignalTime) {
        state.lastSignalTime = time;
        return;
      }

      let delta = time - state.lastSignalTime;
      state.lastSignalTime = time;

      if (delta > 140) {
        delta = state.sampleStepMs;
      }

      state.signalAccumulator += delta;

      while (state.signalAccumulator >= state.sampleStepMs) {
        stepSignals(1);
        state.signalAccumulator -= state.sampleStepMs;
      }
    }

    function drawBufferedTrace(buffer, centerY, color, width) {
      if (!buffer.length) {
        return;
      }

      context.beginPath();
      context.lineWidth = width;
      context.strokeStyle = color;

      buffer.forEach(function(value, index) {
        const x = state.tracks.left + index * state.sampleSpacingPx;
        const y = centerY + value;
        if (index === 0) {
          context.moveTo(x, y);
        } else {
          context.lineTo(x, y);
        }
      });

      context.stroke();
    }

    function drawEcgTrack() {
      drawTrackLabel("Heart beat", state.tracks.ecgY);
      drawGuide(state.tracks.ecgY);
      drawBufferedTrace(state.ecgBuffer, state.tracks.ecgY, "rgba(243, 181, 107, 0.86)", 2.05);
    }

    function drawNeuralTrack() {
      drawTrackLabel("Neural", state.tracks.neuralY);
      drawGuide(state.tracks.neuralY);
      drawBufferedTrace(state.neuralBuffer, state.tracks.neuralY, "rgba(101, 214, 211, 0.82)", 1.7);
    }

    function draw(time) {
      advanceSignals(time);
      context.clearRect(0, 0, state.tracks.width, state.tracks.height);

      const gazePosition = updateGaze(time);
      drawGazeTrack(time, gazePosition);
      drawEcgTrack();
      drawNeuralTrack();

      requestAnimationFrame(draw);
    }

    resize();
    window.addEventListener("resize", resize);
    requestAnimationFrame(draw);
  });
</script>

## Rhythms of Cognition
The human brain accomplishes something we rarely pause to appreciate: it supports a wide range of cognitive functions including perception, attention, memory and speech, within a single, limited neural system. How does the brain coordinate neural computation in time to support these functions? My research suggests one solution is to align neural activity with bodily rhythms, such as eye movements, breathing, and the heartbeat. I have found that each shift in gaze is accompanied by a spatially structured neural response that underlie memory encoding. 

Bodily rhythms carry information about when sensory input is likely to arrive and when internal states are likely to change, providing the brain with a reliable temporal scaffold for organizing computation. Combining human electrophysiology with eye-tracking, I aim to characterize how neural dynamics coordinate with bodily rhythms and how this coordination shapes human cognition.

## Neurotechnology for cognitive health
Cognitive dysfunction in aging and acquired brain injury is an increasingly urgent challenge as the population ages. Building on my research into how brain-body coordination shapes cognition, I pursue two complementary translational strategies: (1) neuromodulation to preserve and restore function by engaging existing neural circuitry and (2) brain-computer interface (BCI) to decode, supplement, or replace neural processing for cognition. Toward these goals, I develop computational methods that dissociate neural patterns underlying cognitive processes to support brain computer interface compensating or replacement of impaired function; and study mechanisms of action of several neural stimulation, including transcranial magnetic stimulation (TMS), transcutaneous auricular vagus nerve stimulation (taVNS), and low-intensity focused ultrasound. These efforts define what, when and how to interface with the brain for cognitive health.

<details class="about-fold">
  <summary>Biography</summary>
  <div class="about-fold__content">
    <p>I am a PhD candidate in Biomedical Engineering at <a href="https://bme.washu.edu">Washington University in St. Louis</a> in the labs Eric Leuthardt and Peter Brunner. I am affiliated with <a href="https://www.neurotechcenter.org/people">National Center for Adaptive Neurotechnologies</a> and <a href="https://neurosurgery.wustl.edu/division-of-neurotechnology-2/">Division of Neurotechnology at Washington University in St. Louis</a>. I am a fellow of <a href="https://sites.wustl.edu/tnnt/">Translational Neuroscience and Neurotechnology Training Program</a>.</p>
    <p>Before joining Washington University in St. Louis, I completed an international dual-degree program in engineering between <a href="http://en.sjtu.edu.cn/">Shanghai Jiao Tong University</a> and CentraleSupélec <a href="https://www.universite-paris-saclay.fr/en">(Paris-Saclay University)</a>, earning a B.Eng. from Shanghai Jiao Tong University and a Diplome d'ingenieur from CentraleSupélec, followed by an M.Eng. in Mechanical Engineering from Shanghai Jiao Tong University in 2022. As part of my dual-degree training, I conducted research at <a href="https://sfrsantelyonest.univ-lyon1.fr/centre51-inserm-u1028-cnrs-umr5292.html">INSERM U1028 Centre de Recherche en Neurosciences de Lyon</a> and was a research fellow at <a href="http://www.l2s.centralesupelec.fr/">L2S: Laboratoire des signaux et systemes</a>. Before beginning my doctoral training, I also served as a visiting scholar in the Department of Neurosurgery at Washington University in St. Louis.</p>
  </div>
</details>

<details class="about-fold">
  <summary>Selected Awards</summary>
  <div class="about-fold__content">
    <ul>
      <li>2025 Biomedical Engineering Teaching Award, Washington University in St. Louis</li>
      <li>2024 Best Scientific Abstract Award, Neurocritical Care Society</li>
      <li>2021 China National Scholarship</li>
    </ul>
  </div>
</details>
