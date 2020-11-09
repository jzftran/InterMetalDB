from django import template

register = template.Library()

@register.filter('get_value_from_dict')
def get_value_from_dict(dict_data, key):
    """
    usage example {{ your_dict|get_value_from_dict:your_key }}
    """
    if key:
        return dict_data.get(key)

@register.filter('capitalize')
def capitalize(string):
    output = ''.join(c for c in string if not c.isnumeric())
    return output.capitalize()
