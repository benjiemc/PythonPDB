# CHANGELOG



## v0.1.3 (2024-03-19)

### Fix

* fix: Add specific version to action ([`05b44e9`](https://github.com/benjiemc/PythonPDB/commit/05b44e9aa3c8dab6cf72aca3b0e38a4326b4b28f))

### Unknown

* Merge pull request #32 from benjiemc/hot-fix

fix: Add specific version to action ([`50dc19f`](https://github.com/benjiemc/PythonPDB/commit/50dc19fb4c1526af55019a3ae27374568c4f5eaa))


## v0.1.2 (2024-03-19)

### Documentation

* docs: Move script ([`39a8bf8`](https://github.com/benjiemc/PythonPDB/commit/39a8bf895477c33f0458cbe0fa69ff5771b08080))

* docs: Prepend README to documentation home page ([`e6d4fd0`](https://github.com/benjiemc/PythonPDB/commit/e6d4fd0779212bc845ea9d5b68b0f6f8edf9c000))

* docs: Remove trailing colon in heading titles ([`902c6b9`](https://github.com/benjiemc/PythonPDB/commit/902c6b9366cbccd838df5751e1bb4dfadf85d7ad))

* docs: Add example to documentation ([`386ab1c`](https://github.com/benjiemc/PythonPDB/commit/386ab1c2e189ff6e2acbe34a5bf88bff48593153))

* docs: Add title to example notebook ([`16861c6`](https://github.com/benjiemc/PythonPDB/commit/16861c62aed2e794262857e1b6ffca620af2ed81))

* docs: Fix links ([`09b58f8`](https://github.com/benjiemc/PythonPDB/commit/09b58f8269130521fa6dbce8d1a522067078eebd))

* docs: Add information to README ([`9826eed`](https://github.com/benjiemc/PythonPDB/commit/9826eedd82d172d63d50f1d4275a2f453c71e77b))

### Fix

* fix: Add pandonc to documentation CD pipeline ([`7de110e`](https://github.com/benjiemc/PythonPDB/commit/7de110e3d2b1a238366a9895c259164feb148122))

### Unknown

* Merge pull request #31 from benjiemc/hot-fix

fix: Add pandonc to documentation CD pipeline ([`e54eafe`](https://github.com/benjiemc/PythonPDB/commit/e54eafe700613854ed3a0652e6836fd7825eed48))

* Merge pull request #30 from benjiemc/improve-documentation

Improve documentation ([`7d73dae`](https://github.com/benjiemc/PythonPDB/commit/7d73dae853a57f98a3b201a05ddbe296897cb78f))


## v0.1.1 (2024-03-12)

### Fix

* fix: Fix bug with HETATMs at the end of file

PDB files sometimes have the HETATMs separate to the rest of the chain
residues (see 1ac6 as an example). This was leading to unexpected
behaviour where two chains with the same name where being added to the
model as oppose to one chain with all residues. ([`19bfd4c`](https://github.com/benjiemc/PythonPDB/commit/19bfd4ce26538a0a786e6548183fae0705d072ac))

### Unknown

* Merge pull request #29 from benjiemc/hot-fix-chain-handling

fix: Fix bug with HETATMs at the end of file ([`74bde2a`](https://github.com/benjiemc/PythonPDB/commit/74bde2a03abe7762c19250551a1f0c84a9461d44))


## v0.1.0 (2024-03-01)

### Feature

* feat: add molecular weights information ([`8da1ddd`](https://github.com/benjiemc/PythonPDB/commit/8da1ddd689601cf177e08af79a5e27ba6898da65))

### Unknown

* Merge pull request #28 from benjiemc/molecular-weights

feat: add molecular weights information ([`53ea0de`](https://github.com/benjiemc/PythonPDB/commit/53ea0de1385029cf06e79f66a4faec897cb4c650))

* Minor fixes ([`73baf0d`](https://github.com/benjiemc/PythonPDB/commit/73baf0d91c69972b40522535a2b372fe2421af33))

* Merge pull request #27 from benjiemc/blogpig-example

Add pandas example ([`e2be8f3`](https://github.com/benjiemc/PythonPDB/commit/e2be8f324d7b472648b6a411accd3df8ed41f3c5))

* Add pandas example ([`7938e67`](https://github.com/benjiemc/PythonPDB/commit/7938e6722e56e941d0d26ccced38e7e4c0ea7f96))


## v0.0.6 (2023-12-15)

### Fix

* fix: add upload to gh artifacts ([`9d4a237`](https://github.com/benjiemc/PythonPDB/commit/9d4a237add60af0a3747cd7f117940643ecc2b2a))

### Unknown

* Remove unused option ([`6cac1d2`](https://github.com/benjiemc/PythonPDB/commit/6cac1d2cf2e6c3dbb51857ca58d430ff72c7f7de))


## v0.0.5 (2023-12-15)

### Fix

* fix: missing src/ path ([`da83ce8`](https://github.com/benjiemc/PythonPDB/commit/da83ce82a0ab99180fc3e2a206cdc26f97bec6ee))

* fix: semantic versioning variables ([`d25c02a`](https://github.com/benjiemc/PythonPDB/commit/d25c02acd06d1d317b73a7d445fd0e73b38d20e0))


## v0.0.4 (2023-12-15)

### Fix

* fix: add check to workflow step ([`9453a62`](https://github.com/benjiemc/PythonPDB/commit/9453a6222025149fa6e075e7c0d9e96c0e6d95bd))

### Unknown

* Merge pull request #26 from benjiemc/publishing

fix: add check to workflow step ([`2d03409`](https://github.com/benjiemc/PythonPDB/commit/2d03409f2b08587642c51772c35848adb550af9d))

* Merge pull request #25 from benjiemc/publishing

fix: Move build command ([`a7be9d3`](https://github.com/benjiemc/PythonPDB/commit/a7be9d351f465db5bef0ea819a125894f6e29a24))

* Move build command ([`f4bea70`](https://github.com/benjiemc/PythonPDB/commit/f4bea7046d1fd9508d17d64178493f8d6c1870f0))


## v0.0.3 (2023-12-15)

### Fix

* fix: pypi release script ([`80f1093`](https://github.com/benjiemc/PythonPDB/commit/80f1093029809c1162bdea6f1181d7d62272e869))

### Unknown

* Merge pull request #24 from benjiemc/publishing

fix: pypi release script ([`d7ac466`](https://github.com/benjiemc/PythonPDB/commit/d7ac466da1bedd883910da7cfd1cb05792469360))


## v0.0.2 (2023-12-15)

### Fix

* fix: Fix release workflow ([`f6a4227`](https://github.com/benjiemc/PythonPDB/commit/f6a4227fefc250f40fd3027f737edf2c8e64ac29))

* fix: Add build command ([`908bda2`](https://github.com/benjiemc/PythonPDB/commit/908bda2cc1d3969882f52b0b80b54e760f5dd612))

### Unknown

* Merge pull request #23 from benjiemc/publishing

fix: Fix release workflow ([`28041d4`](https://github.com/benjiemc/PythonPDB/commit/28041d49437ac83e69784b546cb6a348ef8217c9))

* Merge pull request #22 from benjiemc/publishing

feat: Add build command ([`d32b743`](https://github.com/benjiemc/PythonPDB/commit/d32b7432c9495ee484f04c235ed6b14270ffea59))

* Merge pull request #21 from benjiemc/license

Add LICENSE file ([`b535a47`](https://github.com/benjiemc/PythonPDB/commit/b535a4756c7b81c1969ce172354ba0ea18cb5818))

* Add LICENSE file

Closes #14 ([`268f29b`](https://github.com/benjiemc/PythonPDB/commit/268f29bd1edfe944d6d875b93dfb5f044de075c6))

* Merge pull request #20 from benjiemc/versioning

Update variables for semantic versioning ([`f52f7d1`](https://github.com/benjiemc/PythonPDB/commit/f52f7d124d7244567f4d0e7e09e68b448efa10d5))

* Update variables for semantic versioning ([`3f45225`](https://github.com/benjiemc/PythonPDB/commit/3f452253f3cdee33c73c71a510e71117ac27b00d))


## v0.0.1 (2023-10-26)

### Fix

* fix: change branch master -&gt; main ([`a4791dc`](https://github.com/benjiemc/PythonPDB/commit/a4791dc59b64ad587f127bfb350d07fa2d478380))

### Unknown

* Merge pull request #19 from benjiemc/versioning

fix: change branch master -&gt; main ([`3c40e70`](https://github.com/benjiemc/PythonPDB/commit/3c40e707676fbe09674dde997ca3e3a9abe6738b))

* Merge pull request #18 from benjiemc/versioning

Closes #12 ([`6361d46`](https://github.com/benjiemc/PythonPDB/commit/6361d46db116cd7d7c049f3fa3851bfd76e9633b))

* Add semantic versioning to project ([`f644088`](https://github.com/benjiemc/PythonPDB/commit/f6440888759b86a04bc6e8f611babcc1a2e090b0))

* Merge pull request #17 from benjiemc/update-alignment-code

Update alignment code ([`cfe30df`](https://github.com/benjiemc/PythonPDB/commit/cfe30dfb5ee985cef2e8a7b97f8923f3ba03f2f4))

* Refactor aligners

- rename `align` to `align_structures`
- separate code into functions
- standardize API used in tests ([`b220134`](https://github.com/benjiemc/PythonPDB/commit/b22013410829bb3498e6d27736a2880b4b7453ce))

* Add align function for dataframes ([`b3fd3ce`](https://github.com/benjiemc/PythonPDB/commit/b3fd3ceda9c0c9b8be209d7ac240bdf919c75808))

* Merge pull request #16 from benjiemc/pandas-function

Add function to convert PDB file into a dataframe ([`ce2dae4`](https://github.com/benjiemc/PythonPDB/commit/ce2dae4116302c7bbeea285424d86a8289448ec5))

* Add function to convert PDB file into a dataframe

This adds functionality so that a pdb file doesn&#39;t need to be converted
to structure object and then a dataframe.

Closes #15 ([`748ee88`](https://github.com/benjiemc/PythonPDB/commit/748ee889620560014baa675b3bc2f118b5c1f4db))

* Single source package version ([`b79658e`](https://github.com/benjiemc/PythonPDB/commit/b79658ed6fedc581742c29c79653cbbbd3559523))

* Merge pull request #11 from benjiemc/minor-clean-up

Minor clean up ([`ab4b72a`](https://github.com/benjiemc/PythonPDB/commit/ab4b72af20f71000b0f04e3e1dca4478f8a70872))

* Use standard amino acid list ([`b87b81b`](https://github.com/benjiemc/PythonPDB/commit/b87b81bec7ddab1fd1320ccc07aea137bbeb1dfb))

* Add workflow for pull requests

Closes #10 ([`caf2db2`](https://github.com/benjiemc/PythonPDB/commit/caf2db2116e2e572165a252e6929eed9771679aa))

* Fix flake8 issues ([`d2f4510`](https://github.com/benjiemc/PythonPDB/commit/d2f451041e8af09d336b5391db3b5fcc5f250a93))

* Merge pull request #9 from fspoendlin/main

Handling of exceptions ([`3e57580`](https://github.com/benjiemc/PythonPDB/commit/3e5758091b2b84f325f802029a5ae1d3a9820d4d))

* merge with remote ([`76afd64`](https://github.com/benjiemc/PythonPDB/commit/76afd640cc811eb6ca967e9c749c1fe071806edf))

* fixed bugs in structure.splitstate method ([`0cb4302`](https://github.com/benjiemc/PythonPDB/commit/0cb4302d5ddc297ce116e2523b5355bdfe35a0a8))

* Merge branch &#39;benjiemc:main&#39; into main ([`cdf1da5`](https://github.com/benjiemc/PythonPDB/commit/cdf1da54b9068cad00d605fba00bfa8b6cdede69))

* Remove &#39;None&#39; states

The split states algorithm added atoms without alt codes in a residue
that already had other atoms with alt codes to a new state since called
&#39;None&#39;. Now, these are copied to all states. ([`6a5c59d`](https://github.com/benjiemc/PythonPDB/commit/6a5c59d3961fe938b2d56dcc038a757e55df5bf2))

* Merge branch &#39;benjiemc:main&#39; into main ([`69a666b`](https://github.com/benjiemc/PythonPDB/commit/69a666b26476ecf53ecd43e1883b97ab79472958))

* bug fixed ([`d40bbd6`](https://github.com/benjiemc/PythonPDB/commit/d40bbd6a9b32fe3664f47d391c2d59645c389450))

* Update documentation ([`99ab368`](https://github.com/benjiemc/PythonPDB/commit/99ab368b913fdda9985d6327e6743d8d38f9891c))

* Fix failing test and fix flake8 issues ([`96844d9`](https://github.com/benjiemc/PythonPDB/commit/96844d98a4b2ac388ba5e1fdf872850138ce2ef9))

* Merge pull request #8 from fspoendlin/main

modified Structure.split_states() function to group alternate residues ([`61b3f03`](https://github.com/benjiemc/PythonPDB/commit/61b3f0340d0be9b03041ba8bc287cde7cb9b8166))

* modified Structure.split_states() function to group alternate residues ([`0b308bc`](https://github.com/benjiemc/PythonPDB/commit/0b308bcf3501bb77d2f28dafcc821ca67bb12645))

* Merge pull request #6 from benjiemc/aligner

Aligner ([`6fc2a77`](https://github.com/benjiemc/PythonPDB/commit/6fc2a77a2fae6d38aa81ec79dbf9eaa46c86fa09))

* Update align function to take coordinates ([`434a47b`](https://github.com/benjiemc/PythonPDB/commit/434a47ba57ba864a29fe8aa655dd5c2160d8c23e))

* Added sequence alignment function ([`90ca7ff`](https://github.com/benjiemc/PythonPDB/commit/90ca7ffe79f973f7d97223d8bb7c87e45d161723))

* Aligner for structures ([`873c5b6`](https://github.com/benjiemc/PythonPDB/commit/873c5b63ae0b095fe28fdcbabb119ae879c30cb2))

* Add rmsd function ([`478a2f2`](https://github.com/benjiemc/PythonPDB/commit/478a2f2a20042f81f9b2781445c437cbf9438212))

* Add get_coordinates method to entities ([`217ba7e`](https://github.com/benjiemc/PythonPDB/commit/217ba7e9dc7bf3288b6af0ffddbc4a44cbef8c84))

* Clean up types and docs ([`6aca9c4`](https://github.com/benjiemc/PythonPDB/commit/6aca9c4a0a53b409166b48db37505972adedefa8))

* Create parent class that entities can inherent from

Adds a parent class `Entity` that pulls out the logic from shared
properties of Atoms, Residues, Chains, Models, and Structures.

Closes #2 ([`3312451`](https://github.com/benjiemc/PythonPDB/commit/331245131754c4b5c666af00194a520fbe3595ac))

* Add parsing layer

The purpose of this layer is to remove the logic of structure building
from the different formats. Now there is code that converts the
different possible formats into a standard format, and code that builds
the structure from this standard format. The API remains unchanged after
this.

Closes #1 ([`d9fcf8f`](https://github.com/benjiemc/PythonPDB/commit/d9fcf8f04a8d97b3d46ec36454f6fd8def415511))

* Fix typo in type in documentation

Closes #4 ([`def5510`](https://github.com/benjiemc/PythonPDB/commit/def551079f3e00ef744598f14db0daac3c000303))

* Add  functionality to entities

Closes #3 ([`7bd9658`](https://github.com/benjiemc/PythonPDB/commit/7bd9658a895379d92665244f1d99a7beb89d6475))

* Add dehydrate functionality

Also implements remove_atom/remove_residue functionality to entity
clases. ([`434fca1`](https://github.com/benjiemc/PythonPDB/commit/434fca135d03df581f559880644ca30586ecb41d))

* Clarify docstring ([`92ae5a9`](https://github.com/benjiemc/PythonPDB/commit/92ae5a9fe287ee6ce8dbf7167a579df8229ecbb5))

* Add `silent` option to parser to suppress warnings ([`320c5ed`](https://github.com/benjiemc/PythonPDB/commit/320c5ed42d124ceb787d091d3f357f2d13d00220))

* Add warning to identify multiple conformations when alt codes are present ([`b682fc0`](https://github.com/benjiemc/PythonPDB/commit/b682fc08400e4cef6fb1efd1438a53e4e3fba79a))

* Fix weird comment message ([`4cdd09e`](https://github.com/benjiemc/PythonPDB/commit/4cdd09e220090e2a52463525105f88477069f50a))

* Remove typing import ([`694aab1`](https://github.com/benjiemc/PythonPDB/commit/694aab1546de470b7107d874c0cc4d161c0dc3c0))

* Remove unecessary None check ([`d94af24`](https://github.com/benjiemc/PythonPDB/commit/d94af24280c179b1d181e3a3df41638b1a912caf))

* Make entities equality more holistic

- Fix typo for atom positions
- Remove dependency on parent equality
- Add children to equality check
- Add equality to Structure and Models ([`7f288a2`](https://github.com/benjiemc/PythonPDB/commit/7f288a230954abbd4d0e639be1784cc88ee9e755))

* Add split_states functionality to structure ([`eafa9d4`](https://github.com/benjiemc/PythonPDB/commit/eafa9d43ce72ac687df632e4741038d565d8d033))

* Fix issue with matching seq_ids on different chains

This commit fixes cases where atoms aren&#39;t added to the correct residue
if the same seq_id is seen back to back from one chain to another. ([`fdbe672`](https://github.com/benjiemc/PythonPDB/commit/fdbe67292e7dff75d36fcc092d6820f74a3ab579))

* Add model/endmdl to str representation of structure ([`483b299`](https://github.com/benjiemc/PythonPDB/commit/483b299336b7255b8d7a21dde219b5b895cc735e))

* Add copy functionality to all entities ([`6662ac9`](https://github.com/benjiemc/PythonPDB/commit/6662ac9281f60f1bfc2d8bb619e13acb5eea5c6d))

* String -&gt; str ([`b8d7417`](https://github.com/benjiemc/PythonPDB/commit/b8d7417d4daf4b1d6098e2fc887fefef8fccac4e))

* Fix mistake in alt_loc documentation ([`1a39a2f`](https://github.com/benjiemc/PythonPDB/commit/1a39a2f607b5d63848eeb78e41a628ab540c508e))

* Add .DS_Store to .gitignore ([`47adc57`](https://github.com/benjiemc/PythonPDB/commit/47adc5739ff5220ec25aa82c924e6d7cfcbfea89))

* Restructure formatting things as a subpackage ([`48275f6`](https://github.com/benjiemc/PythonPDB/commit/48275f658cd44bebcdabbb0ba45ea71d2099fd08))

* Improve formats documentation ([`1296e91`](https://github.com/benjiemc/PythonPDB/commit/1296e91ee38a8974444745aee93c89d246f9d523))

* Add support for HETATM records ([`dd4e5ad`](https://github.com/benjiemc/PythonPDB/commit/dd4e5adeecf7679753c80a12173e2c283b320e17))

* Remove &#39;alternates&#39; in Residue ([`c51418f`](https://github.com/benjiemc/PythonPDB/commit/c51418fdfc9a38b001bc83a34b90aa50b4bd421e))

* Add link to documentation in README.md ([`667712b`](https://github.com/benjiemc/PythonPDB/commit/667712ba2a97984afc1270cd52ead662355e1440))

* Add documentation workflow ([`36d223e`](https://github.com/benjiemc/PythonPDB/commit/36d223e6370b5a5d51f06611c11e70409ac532aa))

* Add python 3.11 to test suite ([`9f9e95f`](https://github.com/benjiemc/PythonPDB/commit/9f9e95f3fbead94763d502809c5514b83c62f1fb))

* Add test workflow ([`115c811`](https://github.com/benjiemc/PythonPDB/commit/115c811d4bd16a1a0bcac6ca9a1ffedc55152ee6))

* Move files from other project ([`c3e9dac`](https://github.com/benjiemc/PythonPDB/commit/c3e9dac3ca89a0d1b1381010211b6a67876e44b2))

* Initial commit ([`e2204a4`](https://github.com/benjiemc/PythonPDB/commit/e2204a49712995d816d5a842b3c57cece57ed71e))
