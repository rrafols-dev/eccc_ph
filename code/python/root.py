import os
ROOT = os.path.join(os.path.expanduser('~'), 'Dropbox', 'Slums', 'JDE').replace('\\', '/') + '/'

# Create output folder structure (no-op if it already exists)
for d in ['outputs_csv/mig_full/2050/',
          'outputs_csv/mig_full/lowerbound/',
          'outputs_csv/mig_full/upperbound/',
          'tabFig/']:
    os.makedirs(ROOT + d, exist_ok=True)
