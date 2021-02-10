import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
	long_description = fh.read()

setuptools.setup(
	name="acdc-nn", 
	version="0.0.5",
	author="Giovanni Birolo",
	author_email="giovanni.birolo@unito.it", #FIXME
	description="A small example package", #FIXME
	long_description=long_description,
	long_description_content_type="text/markdown",
	url="https://github.com/compbiomed-unito/acdc-nn",
	packages=setuptools.find_packages(),
	python_requires='>=3.6',
	install_requires=[
		'numpy>=1.19.5',
		'pandas>=1.1.5',
		'Biopython>=1.78',
		'tensorflow>=2.3.1',
		'silence_tensorflow>=1.1.1',
	],
	package_data={
		'acdc_nn': ['weights/*']
		#'tests': ['*'],
		#'': ['weights/*', 'tests/*'],
	},
	classifiers=[
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	],
	entry_points={
		'console_scripts': [
			'acdc-nn=acdc_nn.cmd:main',
		],
	}
)

