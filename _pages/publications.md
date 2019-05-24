---
layout: archive
title: "Publications"
permalink: /publications/
author_profile: true
---
For those articles which are succesfully published, you can find them on <u><a href="https://scholar.google.fr/">Google Scholar</a>.</u>
You can also find Gansheng's articles on <a href="https://www.researchgate.net/profile/Gansheng_Tan">Researchgate.</a><br>

Here is only  the article that are under refining.
===============================================


{% comment %}
{% if author.googlescholar %}
  You can also find my articles on <u><a href="{{author.googlescholar}}">my Google Scholar profile</a>.</u>
{% endif %}

{% include base_path %}
{% endcomment %}

{% for post in site.publications reversed %}
  {% include archive-single.html %}
{% endfor %}
