
from raichu.drawing.colours import OUTLINE_COLOURS, FILL_COLOURS, DOMAIN_ABBREVIATIONS


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
    current_x = 30.0
    font_domains = {'family': 'verdana', 'size': (13)}
    circles = []
    texts = []
    cp_positions = []

    for module in cluster.modules:
        for i, domain in enumerate(module.domains):
            abbreviation = DOMAIN_ABBREVIATIONS.get(domain.type.name)
            if abbreviation is None:
                abbreviation = DOMAIN_ABBREVIATIONS.get(domain.domain_name)
                if abbreviation is None:
                    abbreviation = ''

            if domain.supertype.name == 'UNKNOWN' or domain.supertype.name == 'TAILORING':
                current_y = 23.0
            else:
                current_y = 30.0

            circles.append(make_circle(current_x, current_y, domain))

            opacity = 1.0
            if not domain.used:
                opacity = 0.5

            text = f"""<text 
            x="{current_x}" 
            y="{current_y}" 
            opacity="{opacity}" 
            text-anchor="middle"
            font-family="verdana"
            font-size = "{13}">
            <tspan y="{current_y}" dy="0.35em">{abbreviation}</tspan>
            </text>"""

            texts.append(text)
            if domain.supertype.name == "CARRIER" and domain.used:
                cp_positions.append((current_x, current_y))
            current_x += delta_x

        current_x += min_module_padding

    svg = ''

    for i, circle in enumerate(circles):
        text = texts[i]
        svg += f"""<g id="domain-{i}">"""
        svg += circle
        svg += text
        svg += "</g>"

    return svg, cp_positions, current_x






