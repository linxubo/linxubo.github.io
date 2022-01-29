---
title: "Lin Lab - Publications"
layout: gridlay
excerpt: "Lin Lab -- Publications."
sitemap: false
permalink: /publications/
---


# Publications

## Highlights

(For a full list, see [Google Scholar](https://scholar.google.com/citations?user=VSlcvLoAAAAJ&hl=en) or [ResearchGate](https://www.researchgate.net/profile/Xubo_Lin)).

{% for publi in site.data.publist %}

{% if publi.highlight == 1 %}

<div class="col-sm-12 clearfix">
 <div class="well">
  <img src="{{ site.url }}{{ site.baseurl }}/images/Publications/{{ publi.image }}" class="img-responsive" width="15%" style="float: left" />
  <pubtit>{{ publi.title }}</pubtit>
  <p>{{ publi.description }}</p>
  <p><em>{{ publi.authors }}</em></p>
  <p><strong><a href="{{ publi.link.url }}">{{ publi.link.display }}</a></strong></p>
  <p class="text-danger"><strong> {{ publi.news1 }}</strong></p>
  <p> {{ publi.news2 }}</p>
 </div>
</div>
