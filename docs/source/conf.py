import os
import sys
# Aponta para a pasta raiz onde está o código 'bioinf'
sys.path.insert(0, os.path.abspath('../..'))

# --- Informação ---
project = 'AASB Portfolio'
copyright = '2025, Grupo 06'
author = 'Grupo 06'
release = '1.0'

# --- Extensões ---
extensions = [
    'sphinx.ext.autodoc',      # Lê o teu código
    'sphinx.ext.napoleon',     # Lê docstrings estilo Google
    'sphinx.ext.viewcode',     # Cria links para o código fonte
    'sphinx.ext.coverage',
]

# --- Tema Visual ---
templates_path = ['_templates']
exclude_patterns = []
html_theme = 'sphinx_rtd_theme' # Tema profissional
html_static_path = ['_static']
