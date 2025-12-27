---
title: "MSD Lab - Members"
layout: gridlay
excerpt: "MSD Lab - Members"
sitemap: false
permalink: /members/
---

# Group Members

Jump to [Staffs](#staffs), [Students](#students), [Alumni](#alumni).

## Staffs
{% assign number_printed = 0 %}
{% for member in site.data.team_members %}

{% assign even_odd = number_printed | modulo: 2 %}

{% if even_odd == 0 %}
<div class="row">
{% endif %}

<div class="col-sm-6 clearfix">
  <img src="{{ site.url }}{{ site.baseurl }}/images/members/{{ member.photo }}" class="img-responsive" width="20%" style="float: left" />
  <h4>{{ member.name }}</h4>
  <i>{{ member.info }} <br>Email: <{{ member.email }}></i>
  <ul style="overflow: hidden">

  {% if member.number_educ == 1 %}
  <li> {{ member.education1 }} </li>
  {% endif %}

  {% if member.number_educ == 2 %}
  <li> {{ member.education1 }} </li>
  <li> {{ member.education2 }} </li>
  {% endif %}

  {% if member.number_educ == 3 %}
  <li> {{ member.education1 }} </li>
  <li> {{ member.education2 }} </li>
  <li> {{ member.education3 }} </li>
  {% endif %}

  {% if member.number_educ == 4 %}
  <li> {{ member.education1 }} </li>
  <li> {{ member.education2 }} </li>
  <li> {{ member.education3 }} </li>
  <li> {{ member.education4 }} </li>
  {% endif %}

  {% if member.number_educ == 5 %}
  <li> {{ member.education1 }} </li>
  <li> {{ member.education2 }} </li>
  <li> {{ member.education3 }} </li>
  <li> {{ member.education4 }} </li>
  <li> {{ member.education5 }} </li>
  {% endif %}

  </ul>
</div>

{% assign number_printed = number_printed | plus: 1 %}

{% if even_odd == 1 %}
</div>
{% endif %}

{% endfor %}

{% assign even_odd = number_printed | modulo: 2 %}
{% if even_odd == 1 %}
</div>
{% endif %}




## Students
{% assign number_printed = 0 %}
{% for member in site.data.students %}

{% assign even_odd = number_printed | modulo: 2 %}

{% if even_odd == 0 %}
<div class="row">
{% endif %}

<div class="col-sm-6 clearfix">
  <img src="{{ site.url }}{{ site.baseurl }}/images/members/{{ member.photo }}" class="img-responsive" width="20%" style="float: left" />
  <h4>{{ member.name }}</h4>
  <i>{{ member.info }}<br>Email: <{{ member.email }}></i>
  <ul style="overflow: hidden">

  {% if member.number_educ == 1 %}
  <li> {{ member.education1 }} </li>
  {% endif %}

  {% if member.number_educ == 2 %}
  <li> {{ member.education1 }} </li>
  <li> {{ member.education2 }} </li>
  {% endif %}

  {% if member.number_educ == 3 %}
  <li> {{ member.education1 }} </li>
  <li> {{ member.education2 }} </li>
  <li> {{ member.education3 }} </li>
  {% endif %}

  {% if member.number_educ == 4 %}
  <li> {{ member.education1 }} </li>
  <li> {{ member.education2 }} </li>
  <li> {{ member.education3 }} </li>
  <li> {{ member.education4 }} </li>
  {% endif %}

  </ul>
</div>

{% assign number_printed = number_printed | plus: 1 %}

{% if even_odd == 1 %}
</div>
{% endif %}

{% endfor %}

{% assign even_odd = number_printed | modulo: 2 %}
{% if even_odd == 1 %}
</div>
{% endif %}



## Alumni
1. **Xiu Li**, PhD Student (2021/09~2025/12) & Research Assistant (2020/09~2021/08)
2. **Xiaonian Ji**, MS Student (2023/06~2025/06), then move to DPTechnology (深势科技)
3. **Kaidong Lin**, MS Student (2021/09~2024/06), then move to 北京绿野药业
4. **Xiaoqian Lin**, PhD Student (2020/09~2023/12), then move to Shanghai Jiao Tong University as a postdoc
5. **Shuo Ma**, MS Student (2020/09~2023/06), then move to Smart Logic (思朗科技)
6. **Yifei Zhou**, MS Student (2019/09~2022/06), then move to Aisino (航天信息)
7. **Fujia Tian**, Research Assistant (2018.09-2019.07), then move to City University of Hong Kong
