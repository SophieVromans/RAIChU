
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


def draw_rectangle(gene_name, x, y, width, height=12, text_colour='white'):
    rectangle = f"""<rect x="{x}" y="{y}" width="{width}" height="{height}" fill="black" />"""
    text_x = x + width / 2
    text_y = y + height / 2
    text = f"""<text 
            x="{text_x}" 
            y="{text_y}"
            fill="{text_colour}"
            text-anchor="middle"
            font-family="verdana"
            font-size = "{8}">
            <tspan y="{text_y}" dy="0.35em">{gene_name}</tspan>
            </text>"""

    svg_rect = f"""<g id="gene_banner_{gene_name}">"""
    svg_rect += rectangle
    svg_rect += text
    svg_rect += "</g>"

    return svg_rect


def draw_line(start_x, end_x, y=4):
    line = f"""<line x1="{start_x}" x2="{end_x}" y1="{y}" y2="{y}" stroke-width="2" stroke="black" />"""
    return line


def make_module_text(start_x, end_x, y, module_nr):
    text_x = start_x + (end_x - start_x)/2
    text = f"""<text 
                x="{text_x}" 
                y="{y}"
                fill="'black"
                text-anchor="middle"
                font-family="verdana"
                font-size = "{12}">
                <tspan y="{y}" dy="0.35em">Module {module_nr}</tspan>
                </text>"""
    return text


def draw_bubbles(cluster, widths, delta_x=29, bubble_height=80, min_gene_padding=20, min_module_padding=10):
    x = 30.0
    x_bubbles = 30.0

    cp_positions = []
    cp_positions_bubbles = []
    previous_space_right = 0.0
    previous_cp_position = 0.0

    for i, module in enumerate(cluster.modules):
        current_y = bubble_height
        for j, domain in enumerate(module.domains):
            new_gene = False
            if j < len(module.domains) - 1:
                current_gene = domain.gene
                next_gene = module.domains[j + 1].gene
                if current_gene != next_gene:
                    new_gene = True

            if domain.supertype.name == 'UNKNOWN' or domain.supertype.name == 'TAILORING':

                if current_y == bubble_height - 7:
                    level_change = False
                else:
                    level_change = True

                current_y = bubble_height - 7
            else:
                if current_y == bubble_height:
                    level_change = False
                else:
                    level_change = True
                current_y = bubble_height

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
            if new_gene:
                x_bubbles += min_gene_padding
                x += min_gene_padding

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
    lines = []
    current_gene = None
    gene_start = 0
    gene_end = 0
    gene_rects = []
    genes = []
    module_texts = []

    for i, module in enumerate(cluster.modules):
        start_x = current_x - 15
        current_y = bubble_height
        for j, domain in enumerate(module.domains):
            if domain.gene not in genes:
                genes.append(domain.gene)
            if current_gene and domain.gene != current_gene:
                print("===", domain.type.name)
                width = gene_end - gene_start
                gene_rects.append(draw_rectangle(current_gene, gene_start, bubble_height - 75, width))
                gene_start = current_x - 15
            elif not current_gene:
                gene_start = current_x - 15

            current_gene = domain.gene

            new_gene = False
            if j < len(module.domains) - 1:
                next_gene = module.domains[j + 1].gene
                if current_gene != next_gene:
                    new_gene = True

            abbreviation = DOMAIN_ABBREVIATIONS.get(domain.type.name)
            if abbreviation is None:
                abbreviation = DOMAIN_ABBREVIATIONS.get(domain.domain_name)
                if abbreviation is None:
                    abbreviation = ''

            if domain.supertype.name == 'UNKNOWN' or domain.supertype.name == 'TAILORING':
                if current_y == bubble_height - 7:
                    level_change = False
                else:
                    level_change = True

                current_y = bubble_height - 7
            else:
                if current_y == bubble_height:
                    level_change = False
                else:
                    level_change = True
                current_y = bubble_height

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

            gene_end = current_x + 15
            print(gene_end)

            current_x += delta_x
            if new_gene:

                current_x += min_gene_padding



        # Shift modules based on width of structures

        module_shift = 0.0

        if i != len(cluster.modules) - 1:
            module_shift = module_shifts[i + 1]

        end_x = current_x - 15
        lines.append(draw_line(start_x, end_x, y=bubble_height - 26))
        module_texts.append(make_module_text(start_x, end_x, bubble_height - 34, i))

        current_x += (min_module_padding + module_shift)

    if gene_start and gene_end:
        width = gene_end - gene_start
        gene_rects.append((draw_rectangle(current_gene, gene_start, bubble_height - 75, width)))

    svg = ''

    for i, circle in enumerate(circles):
        text = texts[i]
        svg += f"""<g id="domain-{i}">"""
        svg += circle
        svg += text
        svg += "</g>"

    for line in lines:
        svg += line

    for rect in gene_rects:
        svg += rect

    for module_text in module_texts:
        svg += module_text

    return svg, cp_positions, current_x






