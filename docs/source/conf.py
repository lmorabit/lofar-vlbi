# -*- coding: utf-8 -*-
#
# long_baseline_pipeline documentation build configuration file, created by
# sphinx-quickstart on Mon Nov 11 17:22:41 2019.
#

from datetime import datetime

extensions = []

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

# General information about the project.
project = u'long baseline pipeline'
year = datetime.now().year
copyright = u'%d Leah Morabito + others' % year
author = u'Frits Sweijen, Alexander Drabent and Leah Morabito'

version = u'2.0'
release = u'2.0'
language = None

exclude_patterns = []

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
]

pygments_style = 'sphinx'

html_theme = 'alabaster'
html_sidebars = {
    "**": [
        "about.html",
        "navigation.html",
        "searchbox.html",
    ]
}
    
html_static_path = ['_static']
htmlhelp_basename = 'longbaselinedoc'

html_theme_options = {
    "description": "A calibration pipeline for LOFAR's international stations",
    "github_user": "lmorabit",
    "github_repo": "long_baseline_pipeline",
    "fixed_sidebar": True,
    "github_button" : True,
    #"extra_nav_links": True,
    "show_related": True
}
