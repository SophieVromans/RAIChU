
from raichu.drawing.colours import OUTLINE_COLOURS, FILL_COLOURS, DOMAIN_ABBREVIATIONS, TRANSPARENT_OUTLINE_COLOURS, \
    TRANSPARENT_FILL_COLOURS, TRANSPARENT_TEXT_COLOURS, TEXT_COLOUR


def make_circle(x_coord, y_coord, domain):
    """Easy function to draw circle for the domain visualization. Returns
    matplotlib.patches.circle object of a circle with radius 4

    x_coord: int, x-coordinate of the center of the circle to be drawn
    domain_type: str, PKS domain type (ACP, KS, AT, KR, DH or ER)
    """
    if domain.used:
        colour = FILL_COLOURS[domain.type.name]
        outline_colour = OUTLINE_COLOURS[domain.type.name]
    else:
        colour = TRANSPARENT_FILL_COLOURS[domain.type.name]
        outline_colour = TRANSPARENT_OUTLINE_COLOURS[domain.type.name]

    circle = f"""<circle r="14" cx="{int(x_coord)}" cy="{int(y_coord)}" stroke="{outline_colour}"
       fill="{colour}" stroke-width="1"/>"""
    return circle


def draw_bubbles(cluster, widths, delta_x=29, min_module_padding=10):
    x = 30.0
    x_bubbles = 30.0

    cp_positions = []
    cp_positions_bubbles = []
    previous_space_right = 0.0
    previous_cp_position = 0.0

    for i, module in enumerate(cluster.modules):
        current_y = 30.0
        for j, domain in enumerate(module.domains):

            if domain.supertype.name == 'UNKNOWN' or domain.supertype.name == 'TAILORING':
                if current_y == 23.0:
                    level_change = False
                else:
                    level_change = True

                current_y = 23.0
            else:
                if current_y == 30.0:
                    level_change = False
                else:
                    level_change = True
                current_y = 30.0

            if level_change and j != 0:
                x_bubbles -= 1
                x -= 1


            if domain.supertype.name == "CARRIER" and domain.used:
                cp_positions_bubbles.append(x_bubbles)

                current_space_left, current_space_right = widths[i]

                minimum_x = previous_cp_position + current_space_left + previous_space_right

                x = max([minimum_x, x])
                cp_positions.append(x)
                previous_space_right = current_space_right
                previous_cp_position = x

            x_bubbles += delta_x
            x += delta_x

        x += min_module_padding
        x_bubbles += min_module_padding

    module_shifts = []

    for i, cp_position in enumerate(cp_positions):
        module_shift = cp_positions[i] - cp_positions_bubbles[i]
        for shift in module_shifts:
            module_shift -= shift
        module_shifts.append(module_shift)

    if module_shifts:

        current_x = 30.0 + module_shifts[0]

    else:
        current_x = 30.0

    circles = []
    texts = []
    cp_positions = []

    for i, module in enumerate(cluster.modules):
        current_y = 30.0
        for j, domain in enumerate(module.domains):
            abbreviation = DOMAIN_ABBREVIATIONS.get(domain.type.name)
            if abbreviation is None:
                abbreviation = DOMAIN_ABBREVIATIONS.get(domain.domain_name)
                if abbreviation is None:
                    abbreviation = ''

            if domain.supertype.name == 'UNKNOWN' or domain.supertype.name == 'TAILORING':
                if current_y == 23.0:
                    level_change = False
                else:
                    level_change = True

                current_y = 23.0
            else:
                if current_y == 30.0:
                    level_change = False
                else:
                    level_change = True
                current_y = 30.0

            if level_change and j != 0:
                current_x -= 1

            circles.append(make_circle(current_x, current_y, domain))

            text_colour = TEXT_COLOUR
            if not domain.used:
                text_colour = TRANSPARENT_TEXT_COLOURS[domain.type.name]

            text = f"""<text 
            x="{current_x}" 
            y="{current_y}"
            fill="{text_colour}"
            text-anchor="middle"
            font-family="verdana"
            font-size = "{13}">
            <tspan y="{current_y}" dy="0.35em">{abbreviation}</tspan>
            </text>"""

            texts.append(text)
            if domain.supertype.name == "CARRIER" and domain.used:
                cp_positions.append((current_x, current_y))
            current_x += delta_x

        module_shift = 0.0
        if i != len(cluster.modules) - 1:
            module_shift = module_shifts[i + 1]

        current_x += (min_module_padding + module_shift)

    svg = ''

    for i, circle in enumerate(circles):
        text = texts[i]
        svg += f"""<g id="domain-{i}">"""
        svg += circle
        svg += text
        svg += "</g>"

    return svg, cp_positions, current_x






