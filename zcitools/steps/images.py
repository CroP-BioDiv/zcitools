import os.path
from collections import defaultdict
from .step import Step
from ..utils.show import print_table
from ..utils.exceptions import ZCItoolsValueError

# ToDo: very similar to SequencesStep. Stores some files for some identifiers. Make it general?


class ImagesStep(Step):
    """
Stores list images.
List of image identifier are stored in description.yml.
Each image can be stored in one or more files in different formats.
"""
    _STEP_TYPE = 'images'
    _SUPPORTED_TYPES = ('.jpg', '.jpeg', '.png', '.gif', '.tif', '.tiff', '.ps', '.pdf', '.svg')

    def _init_data(self, type_description):
        self._images = dict()  # seq_ident -> list of files
        if type_description:
            existing_images = self._find_existing_images()
            for image_ident in type_description['images']:
                self._images[image_ident] = existing_images.get(image_ident, [])

    def _check_data(self):
        existing_images = self._find_existing_images()
        exist_image_idents = set(existing_images.keys())
        needed_image_idents = set(self._images.keys())
        # Are all sequences presented
        not_exist = needed_image_idents - exist_image_idents
        if not_exist:
            raise ZCItoolsValueError(f"Image data not presented for: {', '.join(sorted(not_exist))}")

        # Is there more sequences
        more_data = exist_image_idents - needed_image_idents
        if more_data:
            raise ZCItoolsValueError(f"Data exists for not listed image(s): {', '.join(sorted(more_data))}")

    def _find_existing_images(self):
        existing_seqs = defaultdict(list)
        for f in self.step_files(not_cached=True):
            ext = os.path.splitext(f)[1]
            if ext in self._SUPPORTED_TYPES:
                existing_seqs[f[:-len(ext)]].append(f)
        return existing_seqs

    # # Set data
    # def add_image_file(self, f):
    #     # Filename is relative inside step directory
    #     seq_ident, ext = os.path.splitext(f)
    #     if ext not in self._SUPPORTED_TYPES:
    #         raise ZCItoolsValueError(f"Extension '{ext}' is not known image format! {f}")

    #     if image_ident in self._images:
    #         # Remove other (old) files of same sequence
    #         for old_f in self._images[image_ident]:
    #             if old_f != f:
    #                 silent_remove_file(self.step_file(old_f))
    #     else:
    #         self._images[image_ident] = [f]
    #     # #
    #     # self.remove_cache_files()

    def set_images(self, idents):
        for iden in idents:
            if iden not in self._images:
                self._images[iden] = []

    # Save/load data
    def save(self, create=True, needs_editing=False):
        # Store description.yml
        self.save_description(dict(images=sorted(self._images)), create=create, needs_editing=needs_editing)

    # Retrieve data methods
    def image_exists(self, ident):
        return bool(self._images.get(ident))

    def all_images(self):
        return self._images.keys()

    # Show data
    def show_data(self, params=None):
        print_table(['image_ident', 'Files'],
                    [[ident, ', '.join(sorted(fs))] for ident, fs in sorted(self._images.items())])
