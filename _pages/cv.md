---
layout: archive
title: "Curriculum Vitae"
permalink: /cv/
author_profile: true
redirect_from:
  - /resume
---

{% include base_path %}

<style>
  .cv-sheet {
    margin-bottom: 2rem;
  }

  .cv-sheet__hero {
    padding: clamp(2rem, 4vw, 3.2rem);
    border-radius: 28px;
    background:
      radial-gradient(circle at top right, rgba(240, 196, 106, 0.2), transparent 28%),
      radial-gradient(circle at bottom left, rgba(59, 130, 246, 0.14), transparent 34%),
      linear-gradient(135deg, #132238 0%, #16314f 55%, #1f6a73 100%);
    color: #f8fafc;
    box-shadow: 0 30px 62px rgba(19, 34, 56, 0.18);
  }

  .cv-sheet__hero h1 {
    margin: 0;
    font-size: clamp(2.4rem, 5vw, 4.2rem);
    line-height: 1;
    color: #fff;
  }

  .cv-sheet__hero p {
    margin: 0.9rem 0 0;
    max-width: 52rem;
    line-height: 1.8;
    color: rgba(248, 250, 252, 0.9);
  }

  .cv-sheet__meta {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 0.9rem 1.2rem;
    margin-top: 1.35rem;
  }

  .cv-sheet__meta-block {
    padding: 0.9rem 1rem;
    border-radius: 16px;
    background: rgba(255, 255, 255, 0.08);
    border: 1px solid rgba(255, 255, 255, 0.12);
  }

  .cv-sheet__meta-label {
    display: block;
    margin-bottom: 0.28rem;
    font-size: 0.75rem;
    font-weight: 700;
    letter-spacing: 0.08em;
    text-transform: uppercase;
    color: rgba(248, 250, 252, 0.72);
  }

  .cv-sheet__actions {
    display: flex;
    flex-wrap: wrap;
    gap: 0.8rem;
    margin-top: 1.3rem;
  }

  .cv-sheet__actions .btn {
    margin: 0;
  }

  .cv-sheet__section {
    margin-top: 1.4rem;
    padding: 1.35rem;
    border-radius: 20px;
    border: 1px solid rgba(15, 23, 42, 0.08);
    background: linear-gradient(180deg, #fffefb 0%, #f5f7fb 100%);
    box-shadow: 0 16px 36px rgba(15, 23, 42, 0.06);
  }

  .cv-sheet__section h2 {
    margin-top: 0;
  }

  .cv-sheet__lede {
    color: #5b6878;
    font-size: 0.95rem;
  }

  .cv-sheet__interactive-shell {
    display: grid;
    grid-template-columns: minmax(260px, 0.92fr) minmax(320px, 1.18fr);
    gap: 1rem;
    margin-top: 1rem;
    align-items: start;
  }

  .cv-sheet__switch-list {
    display: grid;
    gap: 0.85rem;
  }

  .cv-sheet__switch {
    width: 100%;
    padding: 1rem 1.05rem;
    border: 1px solid rgba(15, 23, 42, 0.08);
    border-radius: 16px;
    background: #fff;
    text-align: left;
    cursor: pointer;
    transition: transform 0.18s ease, box-shadow 0.18s ease, border-color 0.18s ease;
  }

  .cv-sheet__switch:hover,
  .cv-sheet__switch:focus,
  .cv-sheet__switch.is-active {
    transform: translateY(-1px);
    border-color: rgba(31, 106, 115, 0.45);
    box-shadow: 0 18px 34px rgba(15, 23, 42, 0.08);
    outline: none;
  }

  .cv-sheet__switch-row {
    display: grid;
    grid-template-columns: minmax(72px, 92px) 1fr;
    gap: 0.9rem;
    align-items: start;
  }

  .cv-sheet__period {
    font-size: 0.9rem;
    font-weight: 700;
    color: #1f6a73;
    letter-spacing: 0.02em;
  }

  .cv-sheet__switch h3 {
    margin: 0;
    font-size: 1rem;
    color: #102033;
  }

  .cv-sheet__switch-subtitle {
    margin-top: 0.28rem;
    color: #556476;
    font-size: 0.93rem;
    line-height: 1.55;
  }

  .cv-sheet__detail {
    align-self: start;
    padding: 1.2rem;
    border-radius: 18px;
    border: 1px solid rgba(31, 106, 115, 0.16);
    background: linear-gradient(180deg, rgba(31, 106, 115, 0.07) 0%, rgba(255, 255, 255, 0.98) 100%);
  }

  .cv-sheet__panel {
    display: none;
  }

  .cv-sheet__panel.is-active {
    display: block;
  }

  .cv-sheet__panel h3 {
    margin: 0;
    font-size: 1.15rem;
  }

  .cv-sheet__panel-meta {
    margin-top: 0.35rem;
    color: #556476;
    font-size: 0.94rem;
    line-height: 1.6;
  }

  .cv-sheet__panel p {
    margin: 0.8rem 0 0;
    line-height: 1.75;
  }

  .cv-sheet__bullets {
    margin: 0.9rem 0 0;
    padding-left: 1.1rem;
  }

  .cv-sheet__bullets li + li {
    margin-top: 0.45rem;
  }

  .cv-sheet__tag-row {
    display: flex;
    flex-wrap: wrap;
    gap: 0.55rem;
    margin-top: 1rem;
  }

  .cv-sheet__tag {
    display: inline-flex;
    align-items: center;
    padding: 0.42rem 0.72rem;
    border-radius: 999px;
    background: #e8f3f4;
    border: 1px solid rgba(31, 106, 115, 0.18);
    color: #0f3b44;
    font-size: 0.86rem;
    font-weight: 600;
    text-decoration: none;
  }

  .cv-sheet__grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
    gap: 1rem;
  }

  .cv-sheet__fold {
    margin-top: 0.95rem;
    border: 1px solid rgba(15, 23, 42, 0.08);
    border-radius: 16px;
    background: #fff;
    overflow: hidden;
  }

  .cv-sheet__fold summary {
    padding: 0.95rem 1.05rem;
    cursor: pointer;
    font-weight: 700;
    list-style: none;
  }

  .cv-sheet__fold summary::-webkit-details-marker {
    display: none;
  }

  .cv-sheet__fold summary::after {
    content: "+";
    float: right;
    color: #1f6a73;
    font-size: 1.1rem;
  }

  .cv-sheet__fold[open] summary::after {
    content: "-";
  }

  .cv-sheet__fold-body {
    padding: 0 1.05rem 1.05rem;
  }

  .cv-sheet__paper-list {
    margin: 0;
    padding-left: 1.2rem;
  }

  .cv-sheet__paper-list li + li {
    margin-top: 0.7rem;
  }

  .cv-sheet__paper-note {
    margin-top: 0.2rem;
    color: #5b6878;
    font-size: 0.92rem;
  }

  .cv-sheet__list {
    margin: 0;
    padding-left: 1.1rem;
  }

  .cv-sheet__list li + li {
    margin-top: 0.55rem;
  }

  .cv-sheet__mini-card {
    padding: 1rem 1.05rem;
    border-radius: 16px;
    background: #fff;
    border: 1px solid rgba(15, 23, 42, 0.08);
  }

  .cv-sheet__mini-card h3 {
    margin: 0 0 0.55rem;
    font-size: 1rem;
  }

  @media (max-width: 56em) {
    .cv-sheet__interactive-shell {
      grid-template-columns: 1fr;
    }
  }
</style>

<div class="cv-sheet" markdown="0">
  <section class="cv-sheet__hero">
    <h1>Gansheng Tan, Ph.D. Candidate</h1>
    <div class="cv-sheet__meta">
      <div class="cv-sheet__meta-block">
        <span class="cv-sheet__meta-label">Date</span>
        March 4, 2026
      </div>
      <div class="cv-sheet__meta-block">
        <span class="cv-sheet__meta-label">Affiliation</span>
        Washington University School of Medicine<br>
        Department of Neurological Surgery
      </div>
      <div class="cv-sheet__meta-block">
        <span class="cv-sheet__meta-label">Address and Email</span>
        660 S Euclid Ave, St. Louis, MO 63110<br>
        <a href="mailto:g.tan@wustl.edu" style="color:#fff;">g.tan@wustl.edu</a>
      </div>
    </div>
    <div class="cv-sheet__actions">
      <a class="btn btn--inverse btn--large" href="{{ base_path }}/research/">Research</a>
      <a class="btn btn--info btn--large" href="https://www.linkedin.com/in/gansheng-tan-390435147/">LinkedIn</a>
      <a class="btn btn--primary btn--large" href="https://scholar.google.com/citations?user=4wBDfs4AAAAJ">Google Scholar</a>
      <a class="btn btn--large" href="https://github.com/GanshengT">GitHub</a>
    </div>
  </section>

  <section class="cv-sheet__section">
    <h2>Education</h2>
    <div class="cv-sheet__interactive-shell">
      <div class="cv-sheet__switch-list">
        <button class="cv-sheet__switch is-active" type="button" data-group="education" data-target="education-phd">
          <div class="cv-sheet__switch-row">
            <div class="cv-sheet__period">2022-</div>
            <div>
              <h3>Ph.D. in Biomedical Engineering</h3>
              <div class="cv-sheet__switch-subtitle">Washington University in St. Louis<br>St. Louis, Missouri, USA</div>
            </div>
          </div>
        </button>
        <button class="cv-sheet__switch" type="button" data-group="education" data-target="education-meng">
          <div class="cv-sheet__switch-row">
            <div class="cv-sheet__period">2019-2022</div>
            <div>
              <h3>M.Eng. in Mechanical Engineering</h3>
              <div class="cv-sheet__switch-subtitle">Shanghai Jiao Tong University<br>Shanghai, China</div>
            </div>
          </div>
        </button>
        <button class="cv-sheet__switch" type="button" data-group="education" data-target="education-double">
          <div class="cv-sheet__switch-row">
            <div class="cv-sheet__period">2015-2019</div>
            <div>
              <h3>Double Degree Program in Engineering</h3>
              <div class="cv-sheet__switch-subtitle">Paris-Saclay University and Shanghai Jiao Tong University<br>France and China</div>
            </div>
          </div>
        </button>
      </div>
      <div class="cv-sheet__detail">
        <article class="cv-sheet__panel is-active" id="education-phd" data-group="education">
          <h3>Washington University in St. Louis</h3>
          <div class="cv-sheet__panel-meta">Ph.D. in Biomedical Engineering | St. Louis, Missouri, USA</div>
          <ul class="cv-sheet__bullets">
            <li>Advisors: Eric C. Leuthardt and Peter Brunner</li>
            <li>cognitive, system, computational, and translational neuroscience</li>
            <li>Thesis: Towards closed-loop neuromodulation for cognition: revealing the role of saccade-related evoked potentials in human memory encoding and elucidating the cognitive and neuroimmune mechanisms of transcutaneous auricular vagus nerve stimulation.</li>
          </ul>
        </article>
        <article class="cv-sheet__panel" id="education-meng" data-group="education">
          <h3>Shanghai Jiao Tong University</h3>
          <div class="cv-sheet__panel-meta">M.Eng. in Mechanical Engineering | Shanghai, China</div>
          <ul class="cv-sheet__bullets">
            <li>Advisor: Honghai Liu</li>
            <li>Graduate training in mechanical, computational, and biomedical engineering.</li>
            <li>Thesis: Cortico-muscular network, a framework for motor rehabilitation after strokeand optimzing transcranial magnetic stimulation</li>
          </ul>
        </article>
        <article class="cv-sheet__panel" id="education-double" data-group="education">
          <h3>Double Degree Program in Engineering</h3>
          <div class="cv-sheet__panel-meta">Diplome d'ingenieur, Paris-Saclay University | B.Eng., Shanghai Jiao Tong University</div>
          <ul class="cv-sheet__bullets">
            <li>Diplome d'ingenieur (postgraduate degree in engineering), Paris-Saclay University, Gif-sur-Yvette, Paris, France.</li>
            <li>B.Eng., Shanghai Jiao Tong University, Shanghai, China.</li>
            <li>Interdisciplinary training in applied mathematics, systems engineering, signal processing, and machine learning.</li>
            <li>Advisor: Antoine Lutz</li>
          </ul>
        </article>
      </div>
    </div>
  </section>

  <section class="cv-sheet__section">
    <h2>Research and Professional Appointments</h2>
    <div class="cv-sheet__interactive-shell">
      <div class="cv-sheet__switch-list">
        <button class="cv-sheet__switch is-active" type="button" data-group="appointments" data-target="appt-wustl-phd">
          <div class="cv-sheet__switch-row">
            <div class="cv-sheet__period">2022-</div>
            <div>
              <h3>Graduate Student Researcher</h3>
              <div class="cv-sheet__switch-subtitle">Department of Biomedical Engineering and Neurosurgery<br>Washington University in St. Louis, Missouri, USA</div>
            </div>
          </div>
        </button>
        <button class="cv-sheet__switch" type="button" data-group="appointments" data-target="appt-wustl-scholar">
          <div class="cv-sheet__switch-row">
            <div class="cv-sheet__period">2021-2022</div>
            <div>
              <h3>Research Scholar</h3>
              <div class="cv-sheet__switch-subtitle">Department of Neurosurgery<br>Washington University School of Medicine, Missouri, USA</div>
            </div>
          </div>
        </button>
        <button class="cv-sheet__switch" type="button" data-group="appointments" data-target="appt-sjtu">
          <div class="cv-sheet__switch-row">
            <div class="cv-sheet__period">2019-2021</div>
            <div>
              <h3>Graduate Student Researcher</h3>
              <div class="cv-sheet__switch-subtitle">Rehabilitation Medicine and Mechanical Engineering<br>Shanghai Jiao Tong University, Shanghai, China</div>
            </div>
          </div>
        </button>
        <button class="cv-sheet__switch" type="button" data-group="appointments" data-target="appt-inserm">
          <div class="cv-sheet__switch-row">
            <div class="cv-sheet__period">2019</div>
            <div>
              <h3>Research Fellow</h3>
              <div class="cv-sheet__switch-subtitle">Lyon Neuroscience Research Center, INSERM<br>Bron, Lyon, France</div>
            </div>
          </div>
        </button>
        <button class="cv-sheet__switch" type="button" data-group="appointments" data-target="appt-l2s">
          <div class="cv-sheet__switch-row">
            <div class="cv-sheet__period">2017-2019</div>
            <div>
              <h3>Graduate Student Researcher</h3>
              <div class="cv-sheet__switch-subtitle">Signals and Systems Laboratory, CNRS<br>Ile-de-France, Paris, France</div>
            </div>
          </div>
        </button>
      </div>
      <div class="cv-sheet__detail">
        <article class="cv-sheet__panel is-active" id="appt-wustl-phd" data-group="appointments">
          <h3>Graduate Student Researcher</h3>
          <div class="cv-sheet__panel-meta">Department of Biomedical Engineering and Neurosurgery </div>
          <ul class="cv-sheet__bullets">
            <li>Approach: multimodal human neurophysiological recording across intracranial and scalp EEG, eye tracking, neuroimaging, single-unit activity, compound action potentials, and ECG in intraoperative, EMU, and ICU settings.</li>
            <li>Bench-to-bedside validation across invasive and non-invasive neuromodulation technologies including VNS, taVNS, intracranial stimulation, temporal interference, and vibrotactile stimulation.</li>
          </ul>
          <div class="cv-sheet__tag-row">
            <a class="cv-sheet__tag" href="{{ base_path }}/publication/2025-11-15-minds-eye-saccade-related-evoked-potentials-support-visual">Saccade-related neural dynamic</a>
            <a class="cv-sheet__tag" href="{{ base_path }}/publication/2025-01-09-the-effect-of-transcutaneous-auricular-vagus-nerve-stimulation">taVNS randomized trial</a>
            <a class="cv-sheet__tag" href="{{ base_path }}/publication/2026-02-01-ultrasound-applications-in-the-treatment-of-major">Ultrasound neuromodulation</a>
            <a class="cv-sheet__tag" href="{{ base_path }}/projects/2025-explainable-machine-learning-cognitive-neuroscience">Explainable ML project</a>
            <a class="cv-sheet__tag" href="{{ base_path }}/projects/2025-Electrode-localization">Electrode localization</a>
          </div>
        </article>
        <article class="cv-sheet__panel" id="appt-wustl-scholar" data-group="appointments">
          <h3>Research Scholar</h3>
          <div class="cv-sheet__panel-meta">Department of Neurosurgery | Advisors: Peter Brunner and Eric C. Leuthardt</div>
          <p>
            Led a prospective clinical trial investigating the anti-inflammatory effects of transcutaneous auricular vagus nerve stimulation in subarachnoid hemorrhage patients and its neurophysiological mechanisms.
          </p>
          <div class="cv-sheet__tag-row">
            <a class="cv-sheet__tag" href="{{ base_path }}/publication/2025-01-09-the-effect-of-transcutaneous-auricular-vagus-nerve-stimulation">taVNS for subarchnoid hemorrhage</a>
            <a class="cv-sheet__tag" href="https://doi.org/10.1371/journal.pone.0301154">NAVSaH protocol</a>
          </div>
        </article>
        <article class="cv-sheet__panel" id="appt-sjtu" data-group="appointments">
          <h3>Graduate Student Researcher</h3>
          <div class="cv-sheet__panel-meta">Department of Rehabilitation Medicine (Ruijin Hospital) and Mechanical Engineering </div>
          <p>
            Pioneered and validated a cortico-muscular network framework for optimizing neuromodulation therapy in motor rehabilitation following stroke.
          </p>
          <div class="cv-sheet__tag-row">
            <a class="cv-sheet__tag" href="{{ base_path }}/publication/2022-04-25-a-framework-for-quantifying-the-effects-of-transcranial">Corticomuscular network</a>
            <a class="cv-sheet__tag" href="{{ base_path }}/publication/2021-06-03-multiscale-transfer-spectral-entropy-for-quantifying-corticomuscular-interaction">Multiscale transfer spectral entropy</a>
          </div>
        </article>
        <article class="cv-sheet__panel" id="appt-inserm" data-group="appointments">
          <h3>Research Fellow</h3>
          <div class="cv-sheet__panel-meta">Lyon Neuroscience Research Center, INSERM | Advisor: Antoine Lutz</div>
          <p>
            Examined how focused attention and open monitoring meditation affect attentional control and perceptual inference using EEG-based mismatch negativity paradigms.
          </p>
          <div class="cv-sheet__tag-row">
            <a class="cv-sheet__tag" href="https://doi.org/10.1109/SMC53654.2022.9945079">EEG meditation classifier paper</a>
          </div>
        </article>
        <article class="cv-sheet__panel" id="appt-l2s" data-group="appointments">
          <h3>Graduate Student Researcher</h3>
          <div class="cv-sheet__panel-meta">Signals and Systems Laboratory, CNRS | Advisor: Antoine Chaillet</div>
          <p>
            Developed an EEG-based neurofeedback system for focused-attention meditation while building the signal-processing and machine-learning foundation for later translational neuroscience work.
          </p>
          <div class="cv-sheet__tag-row">
            <a class="cv-sheet__tag" href="https://doi.org/10.1109/SMC53654.2022.9945079">Focused-attention EEG paper</a>
          </div>
        </article>
      </div>
    </div>
  </section>

  <section class="cv-sheet__section">
    <h2>Publications</h2>
    <p class="cv-sheet__lede"><a href="https://scholar.google.com/citations?user=4wBDfs4AAAAJ">Google Scholar Page</a></p>

    <details class="cv-sheet__fold" open>
      <summary>Peer-reviewed Journal Publications</summary>
      <div class="cv-sheet__fold-body">
        <ol class="cv-sheet__paper-list">
          <li>Tan, G., Huguenard, A. L., Donovan, K. M., Demarest, P., Liu, X., Li, Z., Adamek, M., Lavine, K., Vellimana, A. K., Kummer, T. T., Osbun, J. W., Zipfel, G. J., Brunner, P., &amp; Leuthardt, E. C. (2025). The effect of transcutaneous auricular vagus nerve stimulation on cardiovascular function in subarachnoid hemorrhage patients: A randomized trial. <i>eLife</i>. <a href="https://doi.org/10.7554/elife.100088.3">https://doi.org/10.7554/elife.100088.3</a><div class="cv-sheet__paper-note">Deemed "Important" by eLife editors.</div></li>
          <li>Tan, G., Chen, H., &amp; Leuthardt, E. C. (2025). Ultrasound applications in the treatment of major depressive disorder (MDD): A systematic review of techniques and therapeutic potentials in clinical trials and animal model studies. <i>Neuromodulation: Technology at the Neural Interface</i>. <a href="https://doi.org/10.1016/j.neurom.2025.08.001">https://doi.org/10.1016/j.neurom.2025.08.001</a></li>
          <li>Laurido-Soto, O. J., Tan, G., Nielsen, S. S., Huguenard, A. L., Donovan, K. M., Xu, I., Giles, J., Dhar, R., Adeoye, O., Lee, J.-M., &amp; Leuthardt, E. (2026). Transcutaneous auricular vagus nerve stimulation reduces inflammatory biomarkers after large vessel occlusion stroke: Results of a prospective randomized open-label blinded endpoint trial. <i>Translational Stroke Research</i>, 17(1). <a href="https://doi.org/10.1007/s12975-025-01405-6">https://doi.org/10.1007/s12975-025-01405-6</a></li>
          <li>Donovan, K., Tan, G., Willie, J., Brunner, P., &amp; Leuthardt, E. (2026). Differential gamma responses to transcutaneous auricular vagus nerve stimulation revealed by human intracranial recordings. <i>Neuromodulation: Technology at the Neural Interface</i>. <a href="https://doi.org/10.1016/j.neurom.2025.11.013">https://doi.org/10.1016/j.neurom.2025.11.013</a></li>
          <li>Donovan, K., Adams, J., Park, K., Demarest, P., Tan, G., Willie, J., Brunner, P., Gorlewicz, J., &amp; Leuthardt, E. (2025). Vibrotactile auricular vagus nerve stimulation alters limbic system connectivity in humans: A pilot study. <i>PLOS One</i>, 20(5), e0310917. <a href="https://doi.org/10.1371/journal.pone.0310917">https://doi.org/10.1371/journal.pone.0310917</a></li>
          <li>Tan, G., Adams, J., Donovan, K., Demarest, P., Willie, J. T., Brunner, P., Gorlewicz, J. L., &amp; Leuthardt, E. C. (2024). Does vibrotactile stimulation of the auricular vagus nerve enhance working memory? A behavioral and physiological investigation. <i>Brain Stimulation</i>, 17(2), 460-468. <a href="https://doi.org/10.1016/j.brs.2024.04.002">https://doi.org/10.1016/j.brs.2024.04.002</a></li>
          <li>Huguenard, A., Tan, G., Rivet, D., Gao, F., Johnson, G., Adamek, M., Coxon, A., Kummer, T., Osbun, J., Vellimana, A., Limbrick, D., Zipfel, G., Brunner, P., &amp; Leuthardt, E. (2025). Auricular vagus nerve stimulation for mitigation of inflammation and vasospasm in subarachnoid hemorrhage: a single-institution randomized controlled trial. <i>Journal of Neurosurgery</i>, 1-12. <a href="https://doi.org/10.3171/2024.10.jns241643">https://doi.org/10.3171/2024.10.jns241643</a></li>
          <li>Huguenard, A., Tan, G., Johnson, G., Adamek, M., Coxon, A., Kummer, T., Osbun, J., Vellimana, A., Limbrick Jr, D., Zipfel, G., Brunner, P., &amp; Leuthardt, E. (2024). Non-invasive auricular vagus nerve stimulation for subarachnoid hemorrhage (NAVSaH): protocol for a prospective, triple-blinded, randomized controlled trial. <i>PLOS ONE</i>, 19(8), e0301154. <a href="https://doi.org/10.1371/journal.pone.0301154">https://doi.org/10.1371/journal.pone.0301154</a></li>
          <li>Liu, J., Wang, J., Tan, G., Sheng, Y., Feng, L., Tang, T., Li, X., Xie, Q., Liu, H., &amp; Wei, Y. (2024). A generalized cortico-muscular-cortical network to evaluate the effects of three-week brain stimulation. <i>IEEE Transactions on Biomedical Engineering</i>, 71(1), 195-206. <a href="https://doi.org/10.1109/tbme.2023.3294509">https://doi.org/10.1109/tbme.2023.3294509</a></li>
          <li>Sheng, Y., Wang, J., Tan, G., Chang, H., Xie, Q., &amp; Liu, H. (2024). Muscle synergy plasticity in motor function recovery after stroke. <i>IEEE Transactions on Neural Systems and Rehabilitation Engineering</i>, 32, 1657-1667. <a href="https://doi.org/10.1109/tnsre.2024.3389022">https://doi.org/10.1109/tnsre.2024.3389022</a></li>
          <li>Tan, G., Wang, J., Liu, J., Sheng, Y., Xie, Q., &amp; Liu, H. (2022). A framework for quantifying the effects of transcranial magnetic stimulation on motor recovery from hemiparesis: corticomuscular network. <i>Journal of Neural Engineering</i>, 19(2), 026053. <a href="https://doi.org/10.1088/1741-2552/ac636b">https://doi.org/10.1088/1741-2552/ac636b</a></li>
          <li>Liu, J., Tan, G., Wang, J., Wei, Y., Sheng, Y., Chang, H., Xie, Q., &amp; Liu, H. (2022). Closed-loop construction and analysis of cortico-muscular-cortical functional network after stroke. <i>IEEE Transactions on Medical Imaging</i>, 41(6), 1575-1586. <a href="https://doi.org/10.1109/tmi.2022.3143133">https://doi.org/10.1109/tmi.2022.3143133</a></li>
          <li>Sheng, Y., Tan, G., Liu, J., Chang, H., Wang, J., Xie, Q., &amp; Liu, H. (2022). Upper limb motor function quantification in post-stroke rehabilitation using muscle synergy space model. <i>IEEE Transactions on Biomedical Engineering</i>, 69(10), 3119-3130. <a href="https://doi.org/10.1109/tbme.2022.3161726">https://doi.org/10.1109/tbme.2022.3161726</a></li>
          <li>Tan, G., Xu, K., Liu, J., &amp; Liu, H. (2022). A trend on autism spectrum disorder research: eye tracking-EEG correlative analytics. <i>IEEE Transactions on Cognitive and Developmental Systems</i>, 14(3), 1232-1244. <a href="https://doi.org/10.1109/tcds.2021.3102646">https://doi.org/10.1109/tcds.2021.3102646</a></li>
          <li>Liu, J., Tan, G., Sheng, Y., Wei, Y., &amp; Liu, H. (2022). A novel delay estimation method for improving corticomuscular coherence in continuous synchronization events. <i>IEEE Transactions on Biomedical Engineering</i>, 69(4), 1328-1339. <a href="https://doi.org/10.1109/tbme.2021.3115386">https://doi.org/10.1109/tbme.2021.3115386</a></li>
          <li>Liu, J., Tan, G., Sheng, Y., &amp; Liu, H. (2021). Multiscale transfer spectral entropy for quantifying corticomuscular interaction. <i>IEEE Journal of Biomedical and Health Informatics</i>, 25(6), 2281-2292. <a href="https://doi.org/10.1109/jbhi.2020.3032979">https://doi.org/10.1109/jbhi.2020.3032979</a></li>
          <li>Liu, J., Wang, J., Tan, G., Sheng, Y., Chang, H., Xie, Q., &amp; Liu, H. (2021). Correlation evaluation of functional corticomuscular coupling with abnormal muscle synergy after stroke. <i>IEEE Transactions on Biomedical Engineering</i>, 68(11), 3261-3272. <a href="https://doi.org/10.1109/tbme.2021.3068997">https://doi.org/10.1109/tbme.2021.3068997</a></li>
        </ol>
      </div>
    </details>

    <details class="cv-sheet__fold">
      <summary>Preprints Under Review</summary>
      <div class="cv-sheet__fold-body">
        <ol class="cv-sheet__paper-list">
          <li>Tan, G., Demarest, P., Li, Y., Cho, H., Park, H., Swift, J. R., Inman, C. S., Manns, J. R., Hamann, S. B., Liu, X., Wahlstrom, K. L., Li, Z., Hollearn, M. K., Campbell, J. M., Cettina, P. E., Sivakumar, S. S., Leuthardt, E. C., Willie, J. T., &amp; Brunner, P. (2025). Mind's eye: saccade-related evoked potentials support visual encoding in humans. Cold Spring Harbor Laboratory. <a href="https://doi.org/10.1101/2025.11.11.25339896">https://doi.org/10.1101/2025.11.11.25339896</a></li>
          <li>Tan, G., Laurido-Soto, O., Huguenard, A., Dhar, R., Lee, J., &amp; Leuthardt, E. (2026). Neuroaxonal injury dynamics after acute ischemic stroke and their modulation by transauricular vagus nerve stimulation: results of a prospective randomized open-label blinded endpoint trial NUVISTA.</li>
        </ol>
      </div>
    </details>

    <details class="cv-sheet__fold">
      <summary>Peer-reviewed Conference Papers and Proceedings</summary>
      <div class="cv-sheet__fold-body">
        <ol class="cv-sheet__paper-list">
          <li>Tan, G., Wang, S., Vierge, V., Mu, W., Wang, M., Greco, L., Mounier, H., &amp; Chaillet, A. (2022). An EEG classifier to discriminate between focused attention meditation and problem-solving. <i>2022 IEEE International Conference on Systems, Man, and Cybernetics (SMC)</i>, 1954-1960. <a href="https://doi.org/10.1109/smc53654.2022.9945079">https://doi.org/10.1109/smc53654.2022.9945079</a><div class="cv-sheet__paper-note">Full-length paper.</div></li>
          <li>Tan, G., Wang, J., Liu, J., Sheng, Y., Xie, Q., Brunner, P., &amp; Liu, H. (2022). Towards individualized transcranial magnetic stimulation for motor recovery from hemiparesis: study of corticomuscular network. <i>Neurorehabilitation and Neural Repair</i>, 36(9), NP1-NP38. <a href="https://doi.org/10.1177/15459683221123387">https://doi.org/10.1177/15459683221123387</a></li>
          <li>Huguenard, A., Tan, G., Johnson, G., Adamek, M., Coxon, A., Zipfel, G., Vellimana, A., Brunner, P., &amp; Leuthardt, E. (2023). O-055 Non-invasive auricular vagus nerve stimulation following spontaneous subarachnoid hemorrhage reduces rates of radiographic vasospasm and hospital-acquired infections. <i>SNIS 20th annual meeting oral abstracts</i>, A43-A44. <a href="https://doi.org/10.1136/jnis-2023-snis.55">https://doi.org/10.1136/jnis-2023-snis.55</a></li>
          <li>Liu, J., Tan, G., Sheng, Y., Wang, J., Lu, W., &amp; Liu, H. (2020). Delay estimation for cortical-muscular interaction via the rate of voxel change. <i>2020 IEEE International Conference on Systems, Man, and Cybernetics (SMC)</i>, 3897-3902. <a href="https://doi.org/10.1109/smc42975.2020.9282946">https://doi.org/10.1109/smc42975.2020.9282946</a><div class="cv-sheet__paper-note">Full-length paper.</div></li>
        </ol>
      </div>
    </details>

    <details class="cv-sheet__fold">
      <summary>Patents</summary>
      <div class="cv-sheet__fold-body">
        <ol class="cv-sheet__paper-list">
          <li>Tan, G., Leuthardt, E., Brunner, P. (2025). Brain-Computer Interface for Monitoring and Modulating Neural States Based on Saccade-Related Neural Dynamics. U.S. Provisional Patent Application No. 63/910,588, filed November 3, 2025.</li>
        </ol>
      </div>
    </details>
  </section>

  <section class="cv-sheet__section">
    <h2>Selected Invited Talks and Presentations</h2>
    <details class="cv-sheet__fold" open>
      <summary>Selected Talks</summary>
      <div class="cv-sheet__fold-body">
        <ol class="cv-sheet__paper-list">
          <li>Speaker. 2025 U.S. National Science Foundation Collaborative Research in Computational Neuroscience PI Meeting, UC San Diego. Title: "Saccade-Related Evoked Potentials and Their Role in Human Visual Encoding".</li>
          <li>Invited talk. 2025 Cognitive, Computational, and Systems Neuroscience Pathway Retreat, Washington University in St. Louis. Title: "Understanding and interfacing human cognition: from neurophysiology to neurostimulation".</li>
          <li>Seminar speaker. School of Computer Science and Engineering, Nanyang Technological University, Singapore, 2021. Hosted by Prof. Cuntai Guan. Title: "The interaction between cortical oscillation and muscle synergies".</li>
          <li>Symposium speaker. 30th Annual Graduate Research Symposium, Washington University in St. Louis, 2025. Title: "Transcutaneous Auricular Vagus Nerve Stimulation Reduces Inflammatory Biomarkers and May Improve Outcomes after Large Vessel Occlusion Strokes: Results of the Randomized Clinical Trial NUVISTA".</li>
        </ol>
      </div>
    </details>

    <details class="cv-sheet__fold">
      <summary>Conference Presentations</summary>
      <div class="cv-sheet__fold-body">
        <ol class="cv-sheet__paper-list">
          <li>Presenter, 2025 Society for Neuroscience Annual Meeting, San Diego, USA.
            <div class="cv-sheet__paper-note">Tan, G., et al. The relationship between ongoing oscillatory activity and saccade-related evoked potentials in humans; and co-authored work on face processing in fusiform regions.</div>
          </li>
          <li>Presenter, 2025 American Epilepsy Society Annual Meeting, Atlanta, USA.
            <div class="cv-sheet__paper-note">Co-authored work on distinguishable evoked response profiles in scalp EEG during intraoperative stimulation of the centromedian nucleus of the thalamus.</div>
          </li>
          <li>Presenter, 2024 Society for Neuroscience Annual Meeting, Chicago, USA.
            <div class="cv-sheet__paper-note">Presentations on eye movement and medial temporal lobe activity, taVNS effects during spatial working memory, single-neuron activity during memory retrieval, visual processing localization, cingulate modeling from scalp EEG, and acute effects of cervical VNS.</div>
          </li>
          <li>Presenter, 2023 Society for Neuroscience Annual Meeting, Washington, D.C., USA.
            <div class="cv-sheet__paper-note">Working memory improvement with vibrotactile auricular stimulation and co-authored work on memory-enhancing amygdala stimulation.</div>
          </li>
          <li>Presenter, 2022 Society for Neuroscience Annual Meeting, San Diego, USA.
            <div class="cv-sheet__paper-note">The interaction between cortical oscillation and muscle synergies in patients with hemiparesis.</div>
          </li>
        </ol>
      </div>
    </details>
  </section>

  <section class="cv-sheet__section">
    <h2>Honors and Awards</h2>
    <ul class="cv-sheet__list">
      <li>2026 Translational Neuroscience and Neurotechnology Training Program Fellow, Washington University in St. Louis, United States</li>
      <li>2025 Biomedical Engineering Teaching Award, Washington University in St. Louis, United States</li>
      <li>2024 Best Scientific Abstract Award, Neurocritical Care Society, United States</li>
      <li>2021 China National Scholarship (top 0.5%)</li>
    </ul>
  </section>

  <section class="cv-sheet__grid">
    <article class="cv-sheet__section">
      <h2>Academic Service and Committees</h2>
      <ul class="cv-sheet__list">
        <li>Editor, BMC Biomedical Engineering (2024- )</li>
        <li>Selection Committee Member, IEEE International Conference on Systems, Man, and Cybernetics conference proceedings (2022- )</li>
        <li>Reviewer, American Junior Academy of Science Symposium, American Association for the Advancement of Science</li>
      </ul>
    </article>

    <article class="cv-sheet__section">
      <h2>Reviewer for Scientific Journals</h2>
      <ul class="cv-sheet__list">
        <li>IEEE Transactions on Cybernetics</li>
        <li>IEEE Transactions on Medical Imaging</li>
        <li>IEEE Journal of Biomedical and Health Informatics</li>
        <li>NeuroImage</li>
        <li>Neuromodulation: Technology at the Neural Interface</li>
        <li>iScience</li>
        <li>Progress in Biomedical Engineering</li>
        <li>Journal of Neural Engineering</li>
        <li>Computer Science Review</li>
        <li>NeuroImage: Reports</li>
        <li>PLOS ONE</li>
        <li>Journal of NeuroEngineering and Rehabilitation</li>
        <li>Biomedical Physics and Engineering Express</li>
        <li>Scientific Reports</li>
        <li>Displays</li>
      </ul>
    </article>
  </section>

  <section class="cv-sheet__grid">
    <article class="cv-sheet__section">
      <h2>Open-source Software</h2>
      <ul class="cv-sheet__list">
        <li>Developer of automatic localization of intracranial electrode contacts in human brains. Tan, G. (2024). Intracranial_contact_loc (v1.0). Zenodo. <a href="https://doi.org/10.5281/zenodo.14217838">https://doi.org/10.5281/zenodo.14217838</a></li>
        <li>Developer of automated spike sorting for intracranial recordings. Tan, G. (2025). Human_intracranial_single_unit_sorting (v1.0). Zenodo. <a href="https://doi.org/10.5281/zenodo.15758232">https://doi.org/10.5281/zenodo.15758232</a></li>
        <li>Developer of MNE-Python contributor release record. Larson, E., ... Tan, G., ... (2024). MNE-Python (v1.9.0). Zenodo. <a href="https://doi.org/10.5281/zenodo.14519545">https://doi.org/10.5281/zenodo.14519545</a></li>
        <li>Contributor, Pingouin v0.5.5. Zenodo. <a href="https://doi.org/10.5281/ZENODO.13683424">https://doi.org/10.5281/ZENODO.13683424</a></li>
        <li>Contributor, NeuroKit2. <a href="https://neuropsychology.github.io/NeuroKit/authors.html">https://neuropsychology.github.io/NeuroKit/authors.html</a></li>
      </ul>
    </article>

    <article class="cv-sheet__section">
      <h2>Industry Translation</h2>
      <ul class="cv-sheet__list">
        <li>Transauricular vagus nerve stimulation device for reducing inflammation, Aurenar Inc. (under review for FDA Breakthrough Devices Program).</li>
      </ul>
    </article>
  </section>

  <section class="cv-sheet__section">
    <h2>Teaching, Leadership Activities, and STEM Outreach</h2>
    <ul class="cv-sheet__list">
      <li>2024-2026 Research mentor, Division of Neurotechnology, Washington University School of Medicine. Mentored Tianshu Tan in intracranial electrophysiology; project resulted in a first-author Society for Neuroscience abstract.</li>
      <li>2026 Reviewer, Junior Academy Student Research Colloquium, St. Louis Public Library.</li>
      <li>2025 Teaching team instructor, Young Scientist Program, Washington University School of Medicine. Co-led physiology lessons for K-12 classrooms to broaden STEM participation among historically underserved students.</li>
      <li>2024 Educator, Amazing Brain Carnival, SciFest, St. Louis Science Center. Designed and facilitated interactive neuroscience demonstrations for public audiences.</li>
      <li>2024 Teaching Assistant, BME 5400 Biomedical Data Science, Washington University in St. Louis.</li>
      <li>2021-2022 Instructor, 1SL1000 Convergence, Integration, Probability, and Partial Differential Equation, CentraleSupelec.</li>
      <li>04/2018 - 05/2018 Representative Student, Department of Engineering, University of Cambridge. Exchanged engineering expertise and fostered cultural understanding between the University of Cambridge and CentraleSupelec.</li>
      <li>10/2017 - 06/2018 Project Manager, Tech for Good Initiative, Ile-de-France, France. Led the development of La Condamine, a platform supporting artist-entrepreneurs. <a href="https://lacondamine.org">https://lacondamine.org</a></li>
    </ul>
  </section>

  <section class="cv-sheet__grid">
    <article class="cv-sheet__section">
      <h2>Professional Societies</h2>
      <ul class="cv-sheet__list">
        <li>IEEE</li>
        <li>American Society of Neurorehabilitation</li>
        <li>Society of Neuroscience</li>
      </ul>
    </article>

    <article class="cv-sheet__section">
      <h2>Skills</h2>
      <ul class="cv-sheet__list">
        <li>Software and systems engineering: Python, R, MATLAB, GitHub, Java, C/C++, HTML, CSS</li>
        <li>Biosignal processing and analysis</li>
        <li>Artificial intelligence and advanced statistics</li>
        <li>Clinical and translational research</li>
        <li>Cognitive and systems neuroscience, neurophysiology</li>
        <li>Technological translation and intellectual property management</li>
        <li>Scientific writing and presentation: Adobe Illustrator, MS Office, LaTeX</li>
      </ul>
    </article>

    <article class="cv-sheet__section">
      <h2>Languages</h2>
      <ul class="cv-sheet__list">
        <li>English: fluent</li>
        <li>Chinese: native</li>
        <li>French: proficient</li>
      </ul>
    </article>
  </section>
</div>

<script>
  document.addEventListener("DOMContentLoaded", function() {
    const switches = Array.from(document.querySelectorAll(".cv-sheet__switch"));
    const panels = Array.from(document.querySelectorAll(".cv-sheet__panel"));

    if (!switches.length || !panels.length) return;

    const groups = Array.from(
      new Set(
        switches
          .map(function(button) { return button.dataset.group; })
          .filter(Boolean)
      )
    );

    function showPanel(group, targetId) {
      switches.forEach(function(button) {
        if (button.dataset.group === group) {
          button.classList.toggle("is-active", button.dataset.target === targetId);
        }
      });

      panels.forEach(function(panel) {
        if (panel.dataset.group === group) {
          panel.classList.toggle("is-active", panel.id === targetId);
        }
      });
    }

    groups.forEach(function(group) {
      const groupButtons = switches.filter(function(button) {
        return button.dataset.group === group;
      });

      groupButtons.forEach(function(button) {
        const targetId = button.dataset.target;

        button.addEventListener("mouseenter", function() {
          showPanel(group, targetId);
        });

        button.addEventListener("focus", function() {
          showPanel(group, targetId);
        });

        button.addEventListener("click", function() {
          showPanel(group, targetId);
        });
      });
    });
  });
</script>
