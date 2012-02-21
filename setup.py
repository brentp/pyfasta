from setuptools import setup, find_packages


version = '0.4.5'

setup(name='pyfasta',
      version=version,
      description=\
        "fast, memory-efficient, pythonic (and command-line) access to fasta sequence files",
      url="http://github.com/brentp/pyfasta/",
      long_description=open('README.rst').read() + "\n" + open('CHANGELOG.txt').read(),
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
      keywords='bioinformatics blast fasta',
      author='brentp',
      author_email='bpederse@gmail.com',
      license='MIT',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      package_data={'':['CHANGELOG.txt']},
      include_package_data=True,
      tests_require=['nose'],
      test_suite='nose.collector',
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      scripts=[],
      entry_points={
      'console_scripts': ['pyfasta = pyfasta:main']
      },
  )
