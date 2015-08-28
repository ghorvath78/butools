{%- extends 'basic.tpl' -%}
{% from 'mathjax.tpl' import mathjax %}


{%- block header -%}
<!DOCTYPE html>
<html>
<head>
{%- block html_head -%}
<meta charset="utf-8" />
<title>{{resources['metadata']['name']}}</title>

{% for css in resources.inlining.css -%}
    <style type="text/css">
    {{ css }}
    </style>
{% endfor %}

<!-- Custom stylesheet, it must be in the same directory as the html file -->
<link rel="stylesheet" href="css/bootstrap-theme.min.css">
<link rel="stylesheet" href="css/bootstrap.min.css">

<style type="text/css">
/* Overrides of notebook CSS for static HTML export */
pre {
  font-family: "monospace";
  font-size: 14px;
}
.nbblock {
  padding-right: 5%;
}
</style>


<!-- Loading mathjax macro -->
{{ mathjax() }}
{%- endblock html_head -%}
</head>
{%- endblock header -%}

{% block body %}
<body>
<div class="nbblock">
{{ super() }}
</div>
</body>
{%- endblock body %}

{% block footer %}
</html>
{% endblock footer %}