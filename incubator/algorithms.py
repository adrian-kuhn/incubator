#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import laspy
from skimage.feature import peak_local_max
from skimage.morphology import label, closing, disk, watershed
from skimage.measure import regionprops
from skimage.filters import gaussian
from skimage.restoration import denoise_tv_chambolle
from skimage import io
import rpy2.robjects.packages as rpackages

from . import processing


def install_lidr():
    """
    Setting up r packages to start algorithms in R.
    """
    # import R's utility package
    utils = rpackages.importr('utils')

    # select a mirror for R packages
    utils.chooseCRANmirror(ind=1)  # select the first mirror in the list

    name = "lidR"
    if not rpackages.isinstalled(name):
        utils.install_packages(name)


class WatershedAdaptive(processing.Individual):
    GENE_POOL = {
        "field_denoising_weight": [i for i in numpy.arange(0, 100, 5)],
        "field_sigma": [i for i in numpy.arange(0.0, 10.0, 0.1)],
        "field_truncate": [i for i in numpy.arange(0.0, 30.0, 1.0)],
        "field_min_distance": [i for i in range(1, 30)],
        "field_compactness": [i for i in numpy.arange(20, 300, 10)],
        "forest_denoising_weight": [i for i in numpy.arange(0, 100, 5)],
        "forest_sigma": [i for i in numpy.arange(0.0, 10.0, 0.1)],
        "forest_truncate": [i for i in numpy.arange(0.0, 30.0, 1.0)],
        "forest_min_distance": [i for i in range(1, 30)],
        "forest_compactness": [i for i in numpy.arange(20, 300, 10)],
        "closing_radius": [i for i in numpy.arange(1.0, 5.0, 0.1)],
        "forest_area_threshold": [i for i in numpy.arange(1, 3000, 10)],
    }

    def __init__(self, chromosome=None, assignments=None):
        super().__init__(chromosome, self.GENE_POOL, assignments)

        self.resolution = 4
        self.forest_area_threshold = 0

        self.field_denoising_weight = 0
        self.field_sigma = 0
        self.field_truncate = 0
        self.field_min_distance = 0
        self.field_compactness = 0

        self.forest_denoising_weight = 0
        self.forest_sigma = 0
        self.forest_truncate = 0
        self.forest_min_distance = 0
        self.forest_compactness = 0

        self.closing_radius = 0

    def read_chromosome(self):
        """
        Read the specific chromosome values to class variables.
        """
        self.field_denoising_weight = self.chromosome.genes[0].value
        self.field_sigma = self.chromosome.genes[1].value
        self.field_truncate = self.chromosome.genes[2].value
        self.field_min_distance = self.chromosome.genes[3].value
        self.field_compactness = self.chromosome.genes[4].value
        self.forest_denoising_weight = self.chromosome.genes[5].value
        self.forest_sigma = self.chromosome.genes[6].value
        self.forest_truncate = self.chromosome.genes[7].value
        self.forest_min_distance = self.chromosome.genes[8].value
        self.forest_compactness = self.chromosome.genes[9].value
        self.closing_radius = self.chromosome.genes[10].value
        self.forest_area_threshold = self.chromosome.genes[11].value

    def solve_assignment(self, las_path, las_high_veg_path, dom_path, dtm_path):
        """
        Delineate trees by separating forest and field areas and a watershed segmentation algorithm

        :param las_path: Path to the lidar point cloud stored as las file.
        :type las_path: String
        :param las_high_veg_path: Path to the filtered lidar point cloud stored as las file (high vegetation, class 5)
        :type las_high_veg_path: String
        :param dom_path: Path to the digital surface model of the same area
        :type dom_path: String
        :param dtm_path: Path to the digital terrain model of the same area
        :type dtm_path: String
        :return: Detected tree crowns as a point cloud
        :rtype: :class:`incubator.TreeTops`
        """
        high_veg_las = laspy.file.File(las_high_veg_path, mode="r")

        max_x = numpy.max(high_veg_las.x)
        min_x = numpy.min(high_veg_las.x)
        max_y = numpy.max(high_veg_las.y)
        min_y = numpy.min(high_veg_las.y)
        bin_x = int(round(max_x - min_x))
        bin_y = int(round(max_y - min_y))

        las = laspy.file.File(las_path, mode="r")
        max_x_las = numpy.max(las.x)
        min_x_las = numpy.min(las.x)
        max_y_las = numpy.max(las.y)
        min_y_las = numpy.min(las.y)

        add_left = int(round((min_y - min_y_las) * self.resolution))
        add_top = int(round((min_x - min_x_las) * self.resolution))

        dom = io.imread(dom_path)
        dtm = io.imread(dtm_path)

        # Histogram
        histogram, binx, biny = numpy.histogram2d(high_veg_las.x, high_veg_las.y,
                                                  bins=[(bin_x * self.resolution), (bin_y * self.resolution)])

        # Rotate -90
        histogram = numpy.rot90(histogram)

        # Map the histogram size to the dom
        if bin_x != 100 or bin_y != 100:
            hist = numpy.zeros(dom.shape)
            hist[add_left:add_left + histogram.shape[0], add_top:add_top + histogram.shape[1]] = histogram
        else:
            hist = histogram

        # Closing
        closed = closing(hist, disk(self.closing_radius))

        # Create a mask for regions with trees
        mask = numpy.copy(closed)
        mask[mask != 0] = 1

        veg_dom = numpy.ma.array(dom, mask=(1 - mask).astype(int), fill_value=0).filled()

        # Separating field from forest regions
        regions_field = label(mask)
        regions_forest = numpy.copy(regions_field)
        region_props = regionprops(regions_field, intensity_image=dtm)
        forest_labels = [r.label for r in region_props if
                         r.filled_area / (
                                 self.resolution * self.resolution) > self.forest_area_threshold or r.mean_intensity > 1000]
        regions_forest[numpy.isin(regions_forest, forest_labels, invert=True)] = 0
        regions_field[numpy.isin(regions_field, forest_labels)] = 0

        field = numpy.ma.array(veg_dom, mask=regions_forest, fill_value=0).filled()
        forest = numpy.ma.array(veg_dom, mask=regions_field, fill_value=0).filled()

        trees_field = self.do_watershed(field, self.field_denoising_weight, self.field_sigma, self.field_truncate,
                                        self.field_min_distance, self.field_compactness)
        trees_forest = self.do_watershed(forest, self.forest_denoising_weight, self.forest_sigma, self.forest_truncate,
                                         self.forest_min_distance, self.forest_compactness)

        labels = trees_field + (trees_forest * (numpy.max(trees_field) + 1))

        tree_locations = []
        returned_trees = {}
        for tree in regionprops(labels, intensity_image=dom):
            z = tree.max_intensity

            miny, minx, maxx, maxy = tree.bbox
            loc_y, loc_x = numpy.where(tree.intensity_image == z)
            x = ((minx + loc_x[0]) / self.resolution) + min_x_las
            y = max_y_las - ((miny + loc_y[0]) / self.resolution)

            returned_trees[tree.label] = Tree(tree.label, x, y, z, dom_max=tree.max_intensity)

        for tree in regionprops(labels, intensity_image=dtm):
            t = returned_trees[tree.label]
            ohm = tree.intensity_image[numpy.nonzero(tree.intensity_image)]

            t.dtm_mean = numpy.mean(ohm)
            t.dtm_min = numpy.min(ohm)

            # Correction for min height 3m
            t.height = max(3.0, t.dom_max - t.dtm_mean)
            t.dom_max = max(t.dom_max, t.dtm_mean + t.height)
            t.z = t.dom_max

            tree_locations.append([t.x, t.y, t.z])

        if len(tree_locations) == 0:
            tree_locations.append([0, 0, 0])
        return processing.TreeTops(numpy.vstack(tree_locations))

    def do_watershed(self, img, denoising_weight, sigma, truncate, min_distance, compactness):
        """
        Private function to run the watershed segmentation for a prepared image.

        :param img: The image to run the watershed segmentation
        :type img: Numpy array
        :param denoising_weight: The weight factor for denoising the image before watershed segmentation. See: https://scikit-image.org/docs/dev/api/skimage.restoration.html#skimage.restoration.denoise_tv_chambolle
        :type denoising_weight: Float
        :param sigma: Sigma value for gaussian blurring the image before watershed segmentation. See: https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian
        :type sigma: Float
        :param truncate: Truncation value for gaussian blurring the image before watershed segmentation. See: https://scikit-image.org/docs/dev/api/skimage.filters.html#skimage.filters.gaussian
        :type truncate: Float
        :param min_distance: Minimum distance for local maxia representing tree tops.
        :type min_distance: Float
        :param compactness: Compactness of a watershed basin: See https://scikit-image.org/docs/dev/api/skimage.morphology.html#skimage.morphology.watershed
        :type compactness: Float
        :return: The image containing the labeled areas for individual trees.
        :rtype: Numpy array
        """
        if denoising_weight != 0:
            denoised = denoise_tv_chambolle(img, weight=denoising_weight)
        else:
            denoised = img

        # Gauss filter
        gauss = gaussian(denoised, sigma, mode='constant', cval=0, preserve_range=True, truncate=truncate)

        # Create a mask so the segmentation will only occur on the trees
        mask = numpy.copy(img)
        mask[mask != 0] = 1

        # Local maxima
        local_max = peak_local_max(gauss, indices=False, min_distance=min_distance, exclude_border=False)
        markers = label(local_max)

        # watershed
        labels = watershed(gauss, markers, mask=mask, compactness=compactness)
        return labels

    def mate(self, other):
        """
        Do crossover and mutation with another individual.

        :param other: The other individual to mate this algorithm with
        :type other: :class:`WatershedAdaptive`
        :return: New watershed object build with parameters from parent algorithms
        :rtype: :class:`WatershedAdaptive`
        """
        return WatershedAdaptive(self.chromosome.crossover(other.chromosome), self.assignments)


class Li2012(processing.Individual):
    """
    Delineate Trees according Li et al. (2012).
    Li, W., Guo, Q., Jakubowski, M. K., & Kelly, M. (2012). A new method for segmenting individual trees from the lidar point cloud. Photogrammetric Engineering & Remote Sensing, 78(1), 75-84.
    """
    GENE_POOL = {
        "dt1": [i for i in numpy.arange(0.1, 30.0, 0.1)],
        "dt2": [i for i in numpy.arange(0.1, 30.0, 0.1)],
        "radius": [i for i in numpy.arange(0.0, 10.0, 0.1)],
        "zu": [i for i in numpy.arange(0.1, 40.0, 0.1)],
        "hmin": [i for i in numpy.arange(0.1, 50.0, 0.1)],
        "speed_up": [i for i in range(1, 30)]
    }

    def __init__(self, chromosome=None, assignments=None):
        super().__init__(chromosome, self.GENE_POOL, assignments)

        self.dt1 = 0
        self.dt2 = 0
        self.radius = 0
        self.zu = 0
        self.hmin = 0
        self.speed_up = 0

        install_lidr()

    def read_chromosome(self):
        """
        Read the specific chromosome values to class variables.
        """
        self.dt1 = self.chromosome.genes[0].value
        self.dt2 = self.chromosome.genes[1].value
        self.radius = self.chromosome.genes[2].value
        self.zu = self.chromosome.genes[3].value
        self.hmin = self.chromosome.genes[4].value
        self.speed_up = self.chromosome.genes[5].value

    def solve_assignment(self, las_high_veg_path):
        """
        Delineate trees by starting R environment and using `li2012` function from lidr package.

        :param las_high_veg_path: Path to the filtered lidar point cloud stored as las file.
        :type las_high_veg_path: String
        :return: Detected tree crowns as a point cloud
        :rtype: :class:`incubator.TreeTops`
        """
        lidr = rpackages.importr("lidR")
        stats = rpackages.importr("stats")
        lidr.progress = False
        las = lidr.readLAS(las_high_veg_path, select="xyz", filter="-drop_z_below 0")
        trees = lidr.lastrees(las, lidr.li2012(dt1=self.dt1, dt2=self.dt2, R=self.radius, Zu=self.zu, hmin=self.hmin,
                                               speed_up=self.speed_up))
        trees = stats.na_omit(trees.slots["data"])  # Delete points with no treeID

        tree_data = numpy.transpose(
            numpy.vstack((numpy.array(trees[0]), numpy.array(trees[1]), numpy.array(trees[2]), numpy.array(trees[3]))))
        tree_numbers = numpy.unique(tree_data[:, 3])
        tree_tops = numpy.zeros((len(tree_numbers), 4))
        tree_counter = 0
        for tree_nr in tree_numbers:
            tree = tree_data[numpy.where(tree_data[:, 3] == tree_nr)]
            tree_tops[tree_counter, :] = tree[tree[:, 2].argsort()][-1]
            tree_counter += 1

        tops = tree_tops[:, [0, 1, 2]]

        return processing.TreeTops(tops)

    def mate(self, other):
        """
        Do crossover and mutation with another individual.

        :param other: The other individual to mate this algorithm with
        :type other: :class:`Li2012`
        :return: New watershed object build with parameters from parent algorithms
        :rtype: :class:`Li2012`
        """
        return Li2012(self.chromosome.crossover(other.chromosome), self.assignments)


class Tree:
    """
    Representation of a single tree with multiple parameters.
    """

    def __init__(self, tree_number, x=None, y=None, z=None, dtm_min=None, dtm_mean=None,
                 dom_min=None, dom_mean=None, dom_median=None, dom_max=None, area=None, diameter=None, height=None):
        self.tree_number = tree_number
        self.x = x
        self.y = y
        self.z = z
        self.dtm_min = dtm_min
        self.dtm_mean = dtm_mean
        self.dom_min = dom_min
        self.dom_mean = dom_mean
        self.dom_median = dom_median
        self.dom_max = dom_max
        self.area = area
        self.diameter = diameter
        self.height = height
