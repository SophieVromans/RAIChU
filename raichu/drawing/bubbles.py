
from raichu.drawing.colours import OUTLINE_COLOURS, FILL_COLOURS, DOMAIN_ABBREVIATIONS

# <g color="green">
#     <rect width="50" height="50" fill="currentcolor" />
#     <circle
#       r="25"
#       cx="70"
#       cy="70"
#       stroke="currentcolor"
#       fill="none"
#       stroke-width="5" />
#   </g>
def make_circle(x_coord, y_coord, domain):
    """Easy function to draw circle for the domain visualization. Returns
    matplotlib.patches.circle object of a circle with radius 4

    x_coord: int, x-coordinate of the center of the circle to be drawn
    domain_type: str, PKS domain type (ACP, KS, AT, KR, DH or ER)
    """
    if domain.used:
        alpha = 1.0
    else:
        alpha = 0.5

    circle = f"""<circle r="14" cx="{int(x_coord)}" cy="{int(y_coord)}" stroke="{OUTLINE_COLOURS[domain.type.name]}"
       fill="{FILL_COLOURS[domain.type.name]}" stroke-width="1" opacity="{alpha}"/>"""
    return circle


def draw_bubbles(cluster, delta_x=28, min_module_padding=10):
    current_x = 0.0
    font_domains = {'family': 'verdana', 'size': (13)}
    circles = []
    texts = []

    for module in cluster.modules:
        for i, domain in enumerate(module.domains):
            abbreviation = DOMAIN_ABBREVIATIONS.get(domain.type.name)
            if abbreviation is None:
                abbreviation = DOMAIN_ABBREVIATIONS.get(domain.domain_name)
                if abbreviation is None:
                    abbreviation = ''

            if domain.supertype.name == 'UNKNOWN' or domain.supertype.name == 'TAILORING':
                current_y = -7.0
            else:
                current_y = 0.0

            circles.append(make_circle(current_x, current_y, domain))

            opacity = 1.0
            if not domain.used:
                opacity = 0.5


            texts.append(f"""<text x="{current_x}" y="{current_y}" opacity="{opacity}" text-anchor="middle" dominant-baseline="central" class="small">{abbreviation}</text>""")

            current_x += delta_x

        current_x += min_module_padding

    svg = """<svg viewBox="-20 -20 550 20" xmlns="http://www.w3.org/2000/svg"><style>
    .small {
      font: 13px sans-serif;
    }
    .heavy {
      font: bold 30px sans-serif;
    }

    /* Note that the color of the text is set with the    *
     * fill property, the color property is for HTML only */
    .Rrrrr {
      font: italic 40px serif;
      fill: red;
    }
  </style>"""

    for i, circle in enumerate(circles):
        text = texts[i]
        svg += f"""<g id="domain-{i}">"""
        svg += circle
        svg += text
        svg += "</g>"

    svg += "</svg>"
    return svg






