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
        self.ax = self.fig.add_subplot(111)

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


class _PlotProtein(_Plot):

    def __init__(self, transcript=None, scale=0, colors=None,
                 rename=None, no_domain_labels=False, *args, **kwargs):
        super(_PlotProtein, self).__init__(*args, **kwargs)
        self.scale = scale
        self.transcript = transcript
        self.colors = colors
        self.rename = rename
        self.no_domain_labels = no_domain_labels

    def _scale_protein(self,protein_length):

        # scale the protein

        if self.scale < protein_length:
            self.normalize = protein_length
        else:
            self.normalize = self.scale

        self.offset = 0.05 + (1.0 - float(protein_length)/self.normalize)*0.45
        self.vertical_offset = 0.15

        assert self.normalize >= protein_length, "length normalization should be >= protein length"

    def _draw_domains(self, domains):
        # plot domains

        for domain in domains:

            if domain[2] == '':
                domain_name = str(domain[0])
            else:
                domain_name = str(domain[2])

            color = '#3385ff'
            if domain_name in self.colors:
                color = colors[domain_name]

            if domain_name in self.rename:
                domain_name = rename[domain_name]

            domain_start = (int(domain[3])/float(self.normalize))*0.9 + self.offset
            domain_end = (int(domain[4])/float(self.normalize))*0.9 + self.offset
            domain_center = (domain_end-domain_start)/2. + domain_start

            self.ax.add_patch(
                patches.Rectangle(
                    (
                        domain_start,
                        0.45+self.vertical_offset,
                    ),
                    domain_end-domain_start,
                    0.1,
                    color=color
                )
            )

            if not self.no_domain_labels:

                self.ax.text(
                    domain_center,
                    0.35+self.vertical_offset,
                    domain_name,
                    horizontalalignment='center',
                    fontsize=self.fontsize
                )

    def _draw_protein_length_markers(self, protein_length):
        # plot protein length markers

        self.protein_frame_length = protein_length/float(self.normalize)*0.9

        self.line_end = protein_length/float(self.normalize)*0.9 + self.offset

        self.ax.text(
            0.5,
            0.01+self.vertical_offset,
            "Amino acid position",
            horizontalalignment='center',
            fontsize=self.fontsize
        )

        self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset+self.protein_frame_length
            ),
            (0.3+self.vertical_offset, 0.3+self.vertical_offset),
            color='black'
            )
        )

        # left marker

        self.left_marker_line = self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset
            ),
            (0.25+self.vertical_offset, 0.3+self.vertical_offset),
            color='black'
            )
        )
        self.left_marker_text = self.ax.text(
            self.offset,
            0.15+self.vertical_offset,
            "0",
            horizontalalignment='center',
            fontsize=self.fontsize
        )

        # right marker

        self.right_marker_line = self.ax.add_line(plt.Line2D(
            (
                self.offset+self.protein_frame_length,
                self.offset+self.protein_frame_length
            ),
            (0.25+self.vertical_offset, 0.3+self.vertical_offset),
            color='black'
            )
        )
        self.right_marker_text = self.ax.text(
            self.offset+self.protein_frame_length,
            0.15+self.vertical_offset,
            str(protein_length),
            horizontalalignment='center',
            fontsize=self.fontsize
        )

    def _draw_main_body(self, name_symbols, name_isoform):
        """
        main protein frame
        """

        self.ax.add_patch(
            patches.Rectangle(
                (self.offset, 0.45+self.vertical_offset),
                self.protein_frame_length,
                0.1,
                fill=False
            )
        )

        self.ax.text(
            0.5,
            0.9,
            name_symbols,
            horizontalalignment='center',
            fontsize=self.fontsize
        )
        self.ax.text(
            0.5,
            0.83,
            name_isoform,
            horizontalalignment='center',
            fontsize=self.fontsize-3
        )


class PlotFusionProtein(_PlotProtein):
    def __init__(self,*args, **kwargs):
        super(PlotFusionProtein, self).__init__(*args, **kwargs)

    def _draw_junction(self):
        # add the junction

        self.ax.add_line(plt.Line2D(
            (
                (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset,
                (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
            ),
            (0.4+self.vertical_offset, 0.6+self.vertical_offset),
            color='black'
            )
        )

        # middle marker, loop until it does not overlap with right marker

        overlaps = True
        line_offset = (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
        text_offset = (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
        junction_label_vertical_offset = 0.0

        rr = self.fig.canvas.get_renderer()

        right_marker_text_box = self.right_marker_text.get_window_extent(renderer=rr)
        left_marker_text_box = self.left_marker_text.get_window_extent(renderer=rr)

        while overlaps:

            # middle_marker_line_1/2/3 are to draw angled line

            middle_marker_line_1 = self.ax.add_line(plt.Line2D(
                (
                    (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset,
                    (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
                ),
                (
                    0.275+self.vertical_offset - junction_label_vertical_offset,
                    0.3+self.vertical_offset
                ),
                color='black'
                )
            )

            middle_marker_line_2 = self.ax.add_line(plt.Line2D(
                (
                    line_offset,
                    (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
                ),
                (
                    0.275+self.vertical_offset-junction_label_vertical_offset,
                    0.275+self.vertical_offset-junction_label_vertical_offset
                ),
                color='black'
                )
            )

            middle_marker_line_3 = self.ax.add_line(plt.Line2D(
                (
                    line_offset,
                    line_offset
                ),
                (
                    0.25+self.vertical_offset-junction_label_vertical_offset,
                    0.275+self.vertical_offset-junction_label_vertical_offset
                ),
                color='black'
                )
            )

            middle_marker_text = self.ax.text(
                text_offset,
                0.15+self.vertical_offset-junction_label_vertical_offset,
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

    def draw(self):
        self._scale_protein(self.transcript.protein_length)
        self._draw_domains(self.transcript.domains['fusion'])
        self._draw_protein_length_markers(self.transcript.protein_length)
        self._draw_junction()

        name_symbols = self.transcript.gene5prime.gene.gene_name + ' - ' + \
            self.transcript.gene3prime.gene.gene_name
        name_isoform = self.transcript.transcript1.id + ' - ' + \
            self.transcript.transcript2.id
        self._draw_main_body(name_symbols, name_isoform)

        self.ax.axis('off')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)


class PlotWTProtein(_PlotProtein):
    def __init__(self, ensembl_transcript, *args, **kwargs):
        super(PlotWTProtein, self).__init__(*args, **kwargs)
        self.ensembl_transcript = ensembl_transcript

    def draw(self):
        self._scale_protein(len(self.ensembl_transcript.coding_sequence)/3)
        self._draw_domains(self.transcript.domains[self.ensembl_transcript.id])
        self._draw_protein_length_markers(len(self.ensembl_transcript.coding_sequence)/3)
        self._draw_main_body(
            self.ensembl_transcript.gene.gene_name,
            self.ensembl_transcript.id
        )

        self.ax.axis('off')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
