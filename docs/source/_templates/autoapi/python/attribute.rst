{# Suppress AutoAPI's introspection-derived ``.. py:attribute::`` directive    #}
{# whenever the same attribute is also documented in a NumPy-style              #}
{# ``Attributes`` section of the enclosing class docstring (Napoleon renders   #}
{# that section as the canonical, indexed, cross-referenceable definition).    #}
{#                                                                             #}
{# Heuristic: if the AutoAPI attribute object has no docstring of its own and  #}
{# is not its own page, skip the introspection block.  This covers the common  #}
{# case (``self.x = …`` in ``__init__`` with no inline doc, plus a NumPy       #}
{# ``Attributes`` entry on the class) and prevents duplicate-object-           #}
{# description warnings while keeping all class attributes that *do* carry    #}
{# their own docstring or have their own AutoAPI page.                         #}
{% if obj.display %}
   {% if obj.docstring or is_own_page %}
      {% if is_own_page %}
{{ obj.id }}
{{ "=" * obj.id | length }}

      {% endif %}
.. py:attribute:: {% if is_own_page %}{{ obj.id }}{% else %}{{ obj.name }}{% endif %}
      {% if obj.annotation is not none %}

   :type: {% if obj.annotation %} {{ obj.annotation | tilde_type }}{% endif %}
      {% endif %}
      {% if obj.value is not none %}

         {% if obj.value.splitlines()|count > 1 %}
   :value: Multiline-String

   .. raw:: html

      <details><summary>Show Value</summary>

   .. code-block:: python

      {{ obj.value|indent(width=6,blank=true) }}

   .. raw:: html

      </details>

         {% else %}
   :value: {{ obj.value|truncate(100) }}
         {% endif %}
      {% endif %}

      {% if obj.docstring %}

   {{ obj.docstring|indent(3) }}
      {% endif %}
   {% endif %}
{% endif %}
