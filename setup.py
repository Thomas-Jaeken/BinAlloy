from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="BinAlloy",  # as it would be imported
            # may include packages/namespaces separated by `.`

            # all sources are compiled into a single binary file
            sources=["BinAlloymodule.c", "functions.c"],
            include_dirs=[
                "/Users/thomasjaeken/Desktop/school/MA2/SP2/CompPhys/BinAlloy/"]
        ),
    ]
)
