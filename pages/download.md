---
layout: page
show_meta: false
title: "Download"
header:
   image_fullwidth: "RNA.jpg"
permalink: "/download/"
---

#### Repository

The __PROBer__ GitHub repository is [here](http://github.com/pachterlab/PROBer).

#### Releases

<table class="table">
    <thead>
	<tr>
	    <th style="text-align: left">Version</th>
      	    <th>Date</th>
      	    <th></th>
      	    <th></th>
      	    <th></th>
    	</tr>
    </thead>

    {% for post in site.categories.releases %}
        <tr>
	    <td>Release notes: <a href="{{ site.url }}{{ site.baseurl }}{{ post.url }}">{{ post.version }}</a></td>
     	    <td><span class="entry-date"><time datetime="{{ post.date | date_to_xmlschema }}">{{ post.date | date: "%B %d, %Y" }}</time></span></td>

            <td><a href="https://github.com/pachterlab/PROBer/releases/download/{{ post.version }}/PROBer_mac-{{ post.version }}.tar.gz">Mac</a></td>
	    <td><a href="https://github.com/pachterlab/PROBer/releases/download/{{ post.version }}/PROBer_linux-{{ post.version }}.tar.gz">Linux</a></td>
	    <td><a href="https://github.com/pachterlab/PROBer/archive/{{ post.version }}.tar.gz">Source</a></td>
	</tr>
    {% endfor %}
</table>

#### License

__PROBer__ is distributed under the [GNU General Public License][1].

[1]: {{ site.url }}/license/
