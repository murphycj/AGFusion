import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.pyplot
import json

# this is so I can plot graphics on a headless server

matplotlib.pyplot.ioff()


class _Plot(object):
    def __init__(self, filename='', height=0, width=0, dpi=0, fontsize=12):

        self.filename = filename
        self.width = width
        self.height = height
        self.dpi = dpi
        self.fontsize = fontsize

        self.fig = plt.figure(
            figsize=(self.width, self.height), dpi=self.dpi, frameon=False
        )

    def save(self):

        self.fig.savefig(
            self.filename,
            dpi=self.dpi,
            bbox_inches='tight'
        )

        plt.close(self.fig)
        plt.clf()


class PlotGene:
    def __init__(self):
        pass


class PlotProtein(_Plot):

    def __init__(self, transcript=None, scale=0, colors=None, rename=None,
                 no_domain_labels=False, *args, **kwargs):
        super(PlotProtein, self).__init__(*args, **kwargs)
        self.scale = scale
        self.transcript = transcript
        self.colors = colors
        self.rename = rename
        self.no_domain_labels = no_domain_labels

    def draw(self, transcript=None):
        """
        Draw an individual figure
        """

        ax = self.fig.add_subplot(111)

        # scale the protein

        if self.scale < self.transcript.protein_length:
            normalize = self.transcript.protein_length
        else:
            normalize = self.scale

        offset = 0.05 + (1.0 - float(self.transcript.protein_length)/normalize)*0.45
        vertical_offset = 0.15

        assert normalize >= self.transcript.protein_length, "length normalization should be >= protein length"

        # plot domains

        for domain in self.transcript.domains:

            if domain[2] == '':
                domain_name = str(domain[0])
            else:
                domain_name = str(domain[2])

            color = '#3385ff'
            if domain_name in self.colors:
                color = colors[domain_name]

            if domain_name in self.rename:
                domain_name = rename[domain_name]

            domain_start = (int(domain[3])/float(normalize))*0.9 + offset
            domain_end = (int(domain[4])/float(normalize))*0.9 + offset
            domain_center = (domain_end-domain_start)/2. + domain_start

            ax.add_patch(
                patches.Rectangle(
                    (
                        domain_start,
                        0.45+vertical_offset,
                    ),
                    domain_end-domain_start,
                    0.1,
                    color=color
                )
            )

            if not self.no_domain_labels:

                ax.text(
                    domain_center,
                    0.35+vertical_offset,
                    domain_name,
                    horizontalalignment='center',
                    fontsize=self.fontsize
                )

        # add the junction

        ax.add_line(plt.Line2D(
            (
                (self.transcript.transcript_protein_junction_5prime/float(normalize))*0.9 + offset,
                (self.transcript.transcript_protein_junction_5prime/float(normalize))*0.9 + offset
            ),
            (0.4+vertical_offset, 0.6+vertical_offset),
            color='black'
            )
        )

        # plot protein length markers

        protein_frame_length = self.transcript.protein_length/float(normalize)*0.9

        line_end = self.transcript.protein_length/float(normalize)*0.9 + offset

        ax.text(
            0.5,
            0.01+vertical_offset,
            "Amino acid position",
            horizontalalignment='center',
            fontsize=self.fontsize
        )

        ax.add_line(plt.Line2D(
            (
                offset,
                offset+protein_frame_length
            ),
            (0.3+vertical_offset, 0.3+vertical_offset),
            color='black'
            )
        )

        # left marker

        left_marker_line = ax.add_line(plt.Line2D(
            (
                offset,
                offset
            ),
            (0.25+vertical_offset, 0.3+vertical_offset),
            color='black'
            )
        )
        left_marker_text = ax.text(
            offset,
            0.15+vertical_offset,
            "0",
            horizontalalignment='center',
            fontsize=self.fontsize
        )

        # right marker

        right_marker_line = ax.add_line(plt.Line2D(
            (
                offset+protein_frame_length,
                offset+protein_frame_length
            ),
            (0.25+vertical_offset, 0.3+vertical_offset),
            color='black'
            )
        )
        right_marker_text = ax.text(
            offset+protein_frame_length,
            0.15+vertical_offset,
            str(self.transcript.protein_length),
            horizontalalignment='center',
            fontsize=self.fontsize
        )

        # middle marker, loop until it does not overlap with right marker

        overlaps = True
        line_offset = (self.transcript.transcript_protein_junction_5prime/float(normalize))*0.9 + offset
        text_offset = (self.transcript.transcript_protein_junction_5prime/float(normalize))*0.9 + offset
        junction_label_vertical_offset = 0.0

        rr = self.fig.canvas.get_renderer()

        right_marker_text_box = right_marker_text.get_window_extent(renderer=rr)
        left_marker_text_box = left_marker_text.get_window_extent(renderer=rr)

        while overlaps:

            # middle_marker_line_1/2/3 are to draw angled line

            middle_marker_line_1 = ax.add_line(plt.Line2D(
                (
                    (self.transcript.transcript_protein_junction_5prime/float(normalize))*0.9 + offset,
                    (self.transcript.transcript_protein_junction_5prime/float(normalize))*0.9 + offset
                ),
                (
                    0.275+vertical_offset - junction_label_vertical_offset,
                    0.3+vertical_offset
                ),
                color='black'
                )
            )

            middle_marker_line_2 = ax.add_line(plt.Line2D(
                (
                    line_offset,
                    (self.transcript.transcript_protein_junction_5prime/float(normalize))*0.9 + offset
                ),
                (
                    0.275+vertical_offset-junction_label_vertical_offset,
                    0.275+vertical_offset-junction_label_vertical_offset
                ),
                color='black'
                )
            )

            middle_marker_line_3 = ax.add_line(plt.Line2D(
                (
                    line_offset,
                    line_offset
                ),
                (
                    0.25+vertical_offset-junction_label_vertical_offset,
                    0.275+vertical_offset-junction_label_vertical_offset
                ),
                color='black'
                )
            )

            middle_marker_text = ax.text(
                text_offset,
                0.15+vertical_offset-junction_label_vertical_offset,
                str(self.transcript.transcript_protein_junction_5prime),
                horizontalalignment='center',
                fontsize=self.fontsize
            )

            # detect if text overlaps

            middle_marker_text_box = middle_marker_text.get_window_extent(
                renderer=rr
            )

            # if overlaps then offset the junction text to the left

            if (right_marker_text_box.fully_overlaps(middle_marker_text_box)) and (left_marker_text_box.fully_overlaps(middle_marker_text_box)):
                junction_label_vertical_offset = junction_label_vertical_offset + 0.01

                middle_marker_line_1.remove()
                middle_marker_line_2.remove()
                middle_marker_line_3.remove()
                middle_marker_text.remove()

            elif right_marker_text_box.fully_overlaps(middle_marker_text_box):
                line_offset = line_offset - 0.01
                text_offset = text_offset - 0.01

                middle_marker_line_1.remove()
                middle_marker_line_2.remove()
                middle_marker_line_3.remove()
                middle_marker_text.remove()
            elif left_marker_text_box.fully_overlaps(middle_marker_text_box):
                line_offset = line_offset + 0.01
                text_offset = text_offset + 0.01

                middle_marker_line_1.remove()
                middle_marker_line_2.remove()
                middle_marker_line_3.remove()
                middle_marker_text.remove()
            else:
                overlaps = False

        # main protein frame

        ax.add_patch(
            patches.Rectangle(
                (offset, 0.45+vertical_offset),
                protein_frame_length,
                0.1,
                fill=False
            )
        )

        name_symbols = self.transcript.gene5prime.gene.gene_name + ' - ' + \
            self.transcript.gene3prime.gene.gene_name
        name_isoform = self.transcript.transcript1.id + ' - ' + \
            self.transcript.transcript2.id

        ax.text(
            0.5,
            0.9,
            name_symbols,
            horizontalalignment='center',
            fontsize=self.fontsize
        )
        ax.text(
            0.5,
            0.83,
            name_isoform,
            horizontalalignment='center',
            fontsize=self.fontsize-3
        )

        ax.axis('off')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

        length_normalize = None
