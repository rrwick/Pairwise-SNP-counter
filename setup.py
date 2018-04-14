#!/usr/bin/env python3
import setuptools


# Set package details
package_name = 'pairwise_snp_counter'
package_description = 'Call confident SNPs in pairwise comparison'
package_version = '0.0.1'
author = 'Holt Lab Coding Group'
license = 'GPLv3'


# Call setup
setuptools.setup(
        name=package_name,
        version=package_version,
        license=license,
        test_suite='tests',
        packages=setuptools.find_packages(),
        scripts=['pairwise_snp_counter.py']
        )
