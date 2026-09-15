"""Frozen iteration-4 stream ordering; append only in a new freeze."""
import numpy as np

STREAM_NAMES = ('field', 'outside_lya', 'outside_cmb', 'cmb_noise',
                'quasar_sampling', 'sightline_selection', 'forest_noise',
                'random_catalogue', 'catalogue_split', 'random_templates',
                # iteration 7 (appended; SeedSequence.spawn keys are by index, so the first ten are unchanged):
                'lowz_fields', 'lowz_sampling', 'lowz_randoms')

def seed_streams(seed):
    children = np.random.SeedSequence(int(seed)).spawn(len(STREAM_NAMES))
    return dict(zip(STREAM_NAMES, (np.random.default_rng(s) for s in children)))
