from setuptools import setup

setup(name='formosa',
      version='0.0.2',
      author="Patrick O'Neill",
      author_email="pon2@umbc.edu",
      py_modules=['formosa','maxent_sampling','uniform_sampling', 'formosa_utils', 'evo_sampling'],
      install_requires=['tqdm']
      )
