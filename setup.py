from setuptools import setup, find_packages, Extension

setup(
    name='spkmeansmodule',
    packages=find_packages(),
    ext_modules=[
        Extension(
            'spkmeansmodule',
            ['common/helper_methods.c', 'kmeans/kmeans.c', 'spkmeans.c', 'spkmeansmodule.c']
        )
    ]
)
