from setuptools import setup

setup(name='formosa',
      version='0.0.1',
      author="Patrick O'Neill",
      author_email="pon2@umbc.edu",
      py_modules=['formosa','maxent_sampling','uniform_sampling', 'formosa_utils'],
      install_requires=['tqdm']
      )
