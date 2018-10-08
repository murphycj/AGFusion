import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.pyplot
from itertools import cycle

# this is so I can plot graphics on a headless server

matplotlib.pyplot.ioff()

HORIZONTAL_LEVELS = [1,2,3,4]


class _Plot(object):
    def __init__(self, filename='', height=0, width=0, dpi=0, fontsize=12,
                 scale=0):

        self.filename = filename
        self.scale = scale
        self.width = width
        self.height = height
        self.dpi = dpi
        self.fontsize = fontsize

        self.fig = plt.figure(
            figsize=(self.width, self.height), dpi=self.dpi, frameon=False
        )
        self.ax = self.fig.add_subplot(111)
        self.rr = self.fig.canvas.get_renderer()

    def save(self):

        self.fig.savefig(
            self.filename,
            dpi=self.dpi,
            bbox_inches='tight'
        )

        plt.close(self.fig)
        plt.clf()

    def _scale(self, seq_length):
        """
        scale the sequence (protein or DNA)
        """

        if self.scale is None or self.scale < seq_length:
            self.normalize = seq_length
        else:
            self.normalize = self.scale

        self.offset = 0.05 + (1.0 - float(seq_length)/self.normalize)*0.45

        assert self.normalize >= seq_length, "length normalization should be >= protein length"


class _PlotExons(_Plot):
    def __init__(self, *args, **kwargs):
        super(_PlotExons, self).__init__(*args, **kwargs)
        self.vertical_offset = 0.15

    def _draw_length_markers(self, basepair_length):
        # plot protein length markers

        self.basepair_length = basepair_length/float(self.normalize)*0.9

        self.line_end = basepair_length/float(self.normalize)*0.9 + self.offset

        self.ax.text(
            0.5,
            0.1,
            "Base pair position (kbp)",
            horizontalalignment='center',
            fontsize=self.fontsize
        )

        self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset+self.basepair_length
            ),
            (0.2+self.vertical_offset, 0.2+self.vertical_offset),
            color='black'
            )
        )

        # left marker

        self.left_marker_line = self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset
            ),
            (0.15+self.vertical_offset, 0.2+self.vertical_offset),
            color='black'
            )
        )
        self.left_marker_text = self.ax.text(
            self.offset,
            0.05+self.vertical_offset,
            "0",
            horizontalalignment='center',
            fontsize=self.fontsize
        )

        # draw markers for increments of 1000 base pairs

        for i in range(1, basepair_length+1):
            if (i % 10000) == 0:
                self.left_marker_line = self.ax.add_line(plt.Line2D(
                    (
                        self.offset+(i/float(self.normalize)*0.9),
                        self.offset+(i/float(self.normalize)*0.9)
                    ),
                    (0.175+self.vertical_offset, 0.2+self.vertical_offset),
                    color='black'
                    )
                )

        # right marker

        self.right_marker_line = self.ax.add_line(plt.Line2D(
            (
                self.offset+self.basepair_length,
                self.offset+self.basepair_length
            ),
            (0.15+self.vertical_offset, 0.2+self.vertical_offset),
            color='black'
            )
        )
        self.right_marker_text = self.ax.text(
            self.offset+self.basepair_length,
            0.05+self.vertical_offset,
            str(basepair_length/1000),
            horizontalalignment='center',
            fontsize=self.fontsize
        )


class PlotWTExons(_PlotExons):
    def __init__(self, ensembl_transcript, *args, **kwargs):
        super(PlotWTExons, self).__init__(*args, **kwargs)
        self.ensembl_transcript = ensembl_transcript

    def _draw_main_body(self, name_symbols, name_isoform):
        """
        main protein frame
        """

        length = (self.ensembl_transcript.end-self.ensembl_transcript.start)/float(self.normalize)*0.9

        self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset+length
            ),
            (0.5, 0.5),
            color='black'
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

    def _draw_exons(self):
        for exon in self.ensembl_transcript.exon_intervals:

            if self.ensembl_transcript.strand == '+':
                start = exon[0] - self.ensembl_transcript.start
                end = exon[1] - self.ensembl_transcript.start
            else:

                # this is so the transcription direction is not plotted
                # in reverse for genes on minus strand

                start = -(exon[1] - self.ensembl_transcript.end)
                end = -(exon[0] - self.ensembl_transcript.end)

            exon_start = (int(start)/float(self.normalize))*0.9 + self.offset
            exon_end = (int(end)/float(self.normalize))*0.9 + self.offset
            exon_center = (exon_end-exon_start)/2. + exon_start

            self.ax.add_patch(
                patches.Rectangle(
                    (
                        exon_start,
                        0.45,
                    ),
                    exon_end-exon_start,
                    0.1,
                    color="black"
                )
            )


    def draw(self):
        self._scale(self.ensembl_transcript.end-self.ensembl_transcript.start)
        self._draw_exons()
        self._draw_length_markers(self.ensembl_transcript.end-self.ensembl_transcript.start)
        self._draw_main_body(
            self.ensembl_transcript.gene.gene_name,
            self.ensembl_transcript.id
        )
        self.ax.axis('off')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)


class PlotFusionExons(_PlotExons):
    def __init__(self, transcript, *args, **kwargs):
        super(PlotFusionExons, self).__init__(*args, **kwargs)
        self.transcript = transcript

    def _draw_fusion_junction(self, junction_location):

        junction_location_norm = junction_location/float(self.normalize)*0.9

        self.ax.add_line(plt.Line2D(
            (
                self.offset+junction_location_norm,
                self.offset+junction_location_norm
            ),
            (0.15+self.vertical_offset, 0.2+self.vertical_offset),
            color='black'
            )
        )
        self.right_marker_text = self.ax.text(
            self.offset+junction_location_norm,
            0.05+self.vertical_offset,
            str(junction_location/1000),
            horizontalalignment='center',
            fontsize=self.fontsize
        )

    def _draw_exons(self):
        for exon in self.transcript.gene5prime_exon_intervals:

            if self.transcript.transcript1.strand == '+':
                start = exon[0] - self.transcript.transcript1.start
                end = exon[1] - self.transcript.transcript1.start
            else:

                # this is so the transcription direction is not plotted
                # in reverse for genes on minus strand

                start = -(exon[1] - self.transcript.transcript1.end)
                end = -(exon[0] - self.transcript.transcript1.end)

            exon_start = (int(start)/float(self.normalize))*0.9 + self.offset
            exon_end = (int(end)/float(self.normalize))*0.9 + self.offset
            exon_center = (exon_end-exon_start)/2. + exon_start

            self.ax.add_patch(
                patches.Rectangle(
                    (
                        exon_start,
                        0.45,
                    ),
                    exon_end-exon_start,
                    0.1,
                    color="black"
                )
            )

        if self.transcript.transcript1.strand == '+':
            distance_to_add = self.transcript.gene5prime.junction - \
                self.transcript.transcript1.start
        else:
            distance_to_add = self.transcript.transcript1.end - \
                self.transcript.gene5prime.junction

        for exon in self.transcript.gene3prime_exon_intervals:

            if self.transcript.transcript2.strand == '+':
                start = exon[0] - self.transcript.gene3prime.junction + \
                    distance_to_add
                end = exon[1] - self.transcript.gene3prime.junction + \
                    distance_to_add
            else:

                # this is so the transcription direction is not plotted
                # in reverse for genes on minus strand

                start = (self.transcript.gene3prime.junction - exon[1]) + \
                    distance_to_add
                end = (self.transcript.gene3prime.junction - exon[0]) + \
                    distance_to_add

            exon_start = (int(start)/float(self.normalize))*0.9 + self.offset
            exon_end = (int(end)/float(self.normalize))*0.9 + self.offset
            exon_center = (exon_end-exon_start)/2. + exon_start

            self.ax.add_patch(
                patches.Rectangle(
                    (
                        exon_start,
                        0.45,
                    ),
                    exon_end-exon_start,
                    0.1,
                    color="red"
                )
            )

    def _draw_main_body(self, name_symbols, name_isoform, length):
        """
        main protein frame
        """

        gene5prime_length = 0
        gene3prime_length = 0

        if self.transcript.transcript1.strand == '+':
            gene5prime_length = (self.transcript.gene5prime.junction -
                                 self.transcript.transcript1.start) \
                                 / float(self.normalize)*0.9
        else:
            gene5prime_length = (self.transcript.transcript1.end -
                                 self.transcript.gene5prime.junction) \
                                 / float(self.normalize)*0.9

        if self.transcript.transcript2.strand == '+':
            gene3prime_length = (self.transcript.transcript2.end -
                                 self.transcript.gene3prime.junction) \
                                 / float(self.normalize)*0.9
        else:
            gene3prime_length = (self.transcript.gene3prime.junction -
                                 self.transcript.transcript2.start) \
                                 / float(self.normalize)*0.9

        self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset + gene5prime_length
            ),
            (0.5, 0.5),
            color='black'
            )
        )
        self.ax.add_line(plt.Line2D(
            (
                self.offset + gene5prime_length,
                self.offset + gene5prime_length + gene3prime_length
            ),
            (0.5, 0.5),
            color='red'
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

    def draw(self):

        if self.transcript.transcript1.strand == '+':
            gene5prime_length = self.transcript.gene5prime.junction - \
                self.transcript.transcript1.start
        else:
            gene5prime_length = self.transcript.transcript1.end - \
                self.transcript.gene5prime.junction

        if self.transcript.transcript2.strand == '+':
            gene3prime_length = self.transcript.transcript2.end - \
                self.transcript.gene3prime.junction
        else:
            gene3prime_length = self.transcript.gene3prime.junction - \
                self.transcript.transcript2.start

        self._scale(gene5prime_length+gene3prime_length)
        self._draw_exons()
        self._draw_length_markers(gene5prime_length+gene3prime_length)
        self._draw_fusion_junction(gene5prime_length)
        self._draw_main_body(
            self.transcript.transcript1.gene.gene_name + '-' +
            self.transcript.transcript2.gene.gene_name,
            self.transcript.transcript1.id + '-' +
            self.transcript.transcript2.id,
            gene5prime_length+gene3prime_length
        )
        self.ax.axis('off')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)


class _PlotProtein(_Plot):

    def __init__(self, transcript=None, colors=None,
                 rename=None, no_domain_labels=False,
                 exclude = [], *args, **kwargs):
        super(_PlotProtein, self).__init__(*args, **kwargs)
        self.transcript = transcript
        self.colors = colors
        self.rename = rename
        self.no_domain_labels = no_domain_labels
        self.exclude = exclude
        self.vertical_offset = 0.55

    def _draw_domains(self, domains):
        # plot domains

        domain_labels = {i:[] for i in HORIZONTAL_LEVELS}
        domain_labels_levels = {}
        domain_label_boxes = {i:[] for i in HORIZONTAL_LEVELS}
        lowest_level_plotted = HORIZONTAL_LEVELS[0]

        domains.sort(key=lambda x: x[3])

        domain_count = 0

        for domain in domains:

            # use domain name if available, otherwise use its ID

            if domain[1] is None:
                domain_name = str(domain[0])
            else:
                domain_name = str(domain[1])

            if domain_name in self.exclude:
                continue

            if domain_name in self.rename:
                domain_name = self.rename[domain_name]

            domain_start = (int(domain[3])/float(self.normalize))*0.9 + self.offset
            domain_end = (int(domain[4])/float(self.normalize))*0.9 + self.offset
            domain_center = (domain_end-domain_start)/2. + domain_start

            if not self.no_domain_labels:

                # for each newly plotted domain, loop through all previous
                # plotted domains and calculated the extent of overlap
                # then horizontally adjust the domain label to be plotted
                # closest to the protein domain structure or be
                # on the same level as the label it overlaps with the least

                #domain_stack_level = cycle()
                overlaps = {i:0.0 for i in HORIZONTAL_LEVELS}
                overlaps_all_levels = True
                min_overlap = [float("inf"),HORIZONTAL_LEVELS[0]]
                # plot domain at 1st level

                for level in HORIZONTAL_LEVELS:

                    level_pos = self.vertical_offset - 0.15 - (level-1.0)*0.1

                    tmp_domain_label = self.ax.text(
                        domain_center,
                        level_pos,
                        domain_name,
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=self.fontsize
                    )
                    tmp_domain_label_box = tmp_domain_label.get_window_extent(renderer=self.rr)

                    #check to see if it overlaps with anything

                    if len(domain_label_boxes[level])>0:
                        max_overlap = max([i.x1 - tmp_domain_label_box.x0 for i in domain_label_boxes[level]])

                        if max_overlap > 0.0:
                            overlaps[level] = max_overlap
                            tmp_domain_label.remove()

                            if max_overlap <= min_overlap[0]:
                                min_overlap = [max_overlap, level]
                        else:
                            domain_labels[level].append(tmp_domain_label)
                            domain_label_boxes[level].append(tmp_domain_label_box)
                            overlaps_all_levels = False

                            domain_labels_levels[domain_count] = level

                            if level > lowest_level_plotted:
                                lowest_level_plotted = level

                            break
                    else:
                        domain_labels[level].append(tmp_domain_label)
                        domain_label_boxes[level].append(tmp_domain_label_box)
                        overlaps_all_levels = False

                        domain_labels_levels[domain_count] = level

                        if level > lowest_level_plotted:
                            lowest_level_plotted = level

                        break

                # if the domain label overlaps with something on all levels
                # then plot it on the level where is overlaps the least

                if overlaps_all_levels:

                    level_pos = self.vertical_offset - 0.15 - (min_overlap[1]-1.0)*0.1

                    tmp_domain_label = self.ax.text(
                        domain_center,
                        level_pos,
                        domain_name,
                        horizontalalignment='center',
                        verticalalignment='center',
                        fontsize=self.fontsize
                    )
                    tmp_domain_label_box = tmp_domain_label.get_window_extent(renderer=self.rr)

                    domain_labels[min_overlap[1]].append(tmp_domain_label)
                    domain_label_boxes[min_overlap[1]].append(tmp_domain_label_box)

                    domain_labels_levels[domain_count] = level

                    if min_overlap[1] > lowest_level_plotted:
                        lowest_level_plotted = min_overlap[1]
            domain_count += 1

        # now we know how many levels of domains labels are needed, so
        # remove all levels, make the correction to self.vertical_offset
        # and replot all labels.

        for level, label in list(domain_labels.items()):
            for ll in label:
                ll.remove()

        self.levels_plotted = HORIZONTAL_LEVELS.index(lowest_level_plotted)
        self.vertical_offset += (0.05 * self.levels_plotted)

        domain_count = 0

        for domain in domains:

            if domain[1] is None:
                domain_name = str(domain[0])
            else:
                domain_name = str(domain[1])

            if domain_name in self.exclude:
                continue

            if domain_name in self.rename:
                domain_name = self.rename[domain_name]

            color = '#3385ff'
            if domain_name in self.colors:
                color = self.colors[domain_name]

            domain_start = (int(domain[3])/float(self.normalize))*0.9 + self.offset
            domain_end = (int(domain[4])/float(self.normalize))*0.9 + self.offset
            domain_center = (domain_end-domain_start)/2. + domain_start

            self.ax.add_patch(
                patches.Rectangle(
                    (
                        domain_start,
                        self.vertical_offset,
                    ),
                    domain_end - domain_start,
                    0.1,
                    color=color
                )
            )

            # fetch the level the domain label was determined it was to be
            # plotted on

            level = domain_labels_levels[domain_count]

            level_pos = self.vertical_offset - 0.15 - (level-1.0)*0.1

            tmp_domain_label = self.ax.text(
                domain_center,
                level_pos,
                domain_name,
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=self.fontsize
            )

            domain_count += 1


    def _draw_protein_length_markers(self, protein_length):
        # plot protein length markers

        self.line_end = protein_length/float(self.normalize)*0.9 + self.offset

        self.ax.text(
            0.5,
            self.vertical_offset - (0.5 + self.levels_plotted * 0.1),
            "Amino acid position",
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=self.fontsize
        )

        self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset+self.protein_frame_length
            ),
            (
                self.vertical_offset - (0.35 + self.levels_plotted * 0.05),
                self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
            ),
            color='black'
            )
        )

        # left marker

        self.left_marker_line = self.ax.add_line(plt.Line2D(
            (
                self.offset,
                self.offset
            ),
            (
                self.vertical_offset - (0.38 + self.levels_plotted * 0.05),
                self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
            ),
            color='black'
            )
        )
        self.left_marker_text = self.ax.text(
            self.offset,
            self.vertical_offset - (0.43 + self.levels_plotted * 0.05),
            "0",
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=self.fontsize
        )

        # draw markers for increments of 100 amino acids

        for i in range(1, protein_length+1):
            if (i % 100) == 0:
                self.left_marker_line = self.ax.add_line(plt.Line2D(
                    (
                        self.offset+(i/float(self.normalize)*0.9),
                        self.offset+(i/float(self.normalize)*0.9)
                    ),
                    (
                        self.vertical_offset - (0.38 + self.levels_plotted * 0.05),
                        self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
                    ),
                    color='black'
                    )
                )

        # right marker

        self.right_marker_line = self.ax.add_line(plt.Line2D(
            (
                self.offset+self.protein_frame_length,
                self.offset+self.protein_frame_length
            ),
            (
                self.vertical_offset - (0.38 + self.levels_plotted * 0.05),
                self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
            ),
            color='black'
            )
        )

        self.right_marker_text = self.ax.text(
            self.offset+self.protein_frame_length,
            self.vertical_offset - (0.43 + self.levels_plotted * 0.05),
            str(protein_length),
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=self.fontsize
        )

    def _draw_main_body(self, name_symbols, name_isoform):
        """
        main protein frame
        """

        self.ax.add_patch(
            patches.Rectangle(
                (self.offset, self.vertical_offset),
                self.protein_frame_length,
                0.1,
                fill=False
            )
        )

        self.ax.text(
            0.5,
            0.95,
            name_symbols,
            horizontalalignment='center',
            fontsize=self.fontsize
        )

        self.ax.text(
            0.5,
            0.88,
            name_isoform,
            horizontalalignment='center',
            fontsize=self.fontsize-3
        )


class PlotFusionProtein(_PlotProtein):
    def __init__(self, *args, **kwargs):
        super(PlotFusionProtein, self).__init__(*args, **kwargs)

    def _draw_junction(self):
        # add the junction

        self.ax.add_line(plt.Line2D(
            (
                (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset,
                (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
            ),
            (
                self.vertical_offset - 0.05,
                self.vertical_offset + 0.15
            ),
            color='black'
            )
        )

        # middle marker, loop until it does not overlap with right marker

        overlaps = True
        line_offset = (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
        text_offset = (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
        junction_label_vertical_offset = 0.0

        right_marker_text_box = self.right_marker_text.get_window_extent(renderer=self.rr)
        left_marker_text_box = self.left_marker_text.get_window_extent(renderer=self.rr)

        while overlaps:

            # middle_marker_line_1/2/3 are to draw angled line

            middle_marker_line_1 = self.ax.add_line(plt.Line2D(
                (
                    (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset,
                    (self.transcript.transcript_protein_junction_5prime/float(self.normalize))*0.9 + self.offset
                ),
                (
                    self.vertical_offset - (0.37 + self.levels_plotted * 0.05) - junction_label_vertical_offset,
                    self.vertical_offset - (0.35 + self.levels_plotted * 0.05)
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
                    self.vertical_offset - (0.37 + self.levels_plotted * 0.05) - junction_label_vertical_offset,
                    self.vertical_offset - (0.37 + self.levels_plotted * 0.05) - junction_label_vertical_offset
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
                    self.vertical_offset - (0.4 + self.levels_plotted * 0.05) - junction_label_vertical_offset,
                    self.vertical_offset - (0.37 + self.levels_plotted * 0.05) - junction_label_vertical_offset
                ),
                color='black'
                )
            )

            middle_marker_text = self.ax.text(
                text_offset,
                self.vertical_offset - (0.45 + self.levels_plotted * 0.05) - junction_label_vertical_offset,
                str(self.transcript.transcript_protein_junction_5prime),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=self.fontsize
            )

            # detect if text overlaps

            middle_marker_text_box = middle_marker_text.get_window_extent(
                renderer=self.rr
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
        self._scale(self.transcript.protein_length)
        self.protein_frame_length = self.transcript.protein_length/float(self.normalize)*0.9
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
        self._scale(len(self.ensembl_transcript.coding_sequence)/3)
        self.protein_frame_length = len(self.ensembl_transcript.coding_sequence)/3/float(self.normalize)*0.9
        self._draw_domains(self.transcript.domains[self.ensembl_transcript.id])
        self._draw_protein_length_markers(int(len(self.ensembl_transcript.coding_sequence)/3))
        self._draw_main_body(
            self.ensembl_transcript.gene.gene_name,
            self.ensembl_transcript.id
        )

        self.ax.axis('off')
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
