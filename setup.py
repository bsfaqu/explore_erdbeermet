import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='erdbeermet',
    version='0.0.4',
    author='David Schaller',
    author_email='sdavid@bioinf.uni-leipzig.de',
    description='R matrix simulation and recognition',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/david-schaller/Erdbeermet',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.16.4',
        'scipy>=1.3.0',
        'matplotlib>=3.0',
   ],
)
setuptools.setup(
    name='explore_erdbeermet',
    version='0.0.1',
    author='Bruno Schmidt',
    author_email='bruno@bioinf.uni-leipzig.de',
    description='somewhat interactive R matrix exploration',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/bsfaqu/explore_erdbeermet',
    packages=setuptools.find_packages('src'),
    package_dir={'': 'src'},
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.16.4',
        'scipy>=1.3.0',
        'matplotlib>=3.0',
   ],
)
