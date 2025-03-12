# <img src="https://polybox.ethz.ch/index.php/s/c8uZqrzwk45wpBx/download?path=%2Fpng&files=zonda-high-resolution-logo-transparent.png" width="110" valign="middle" alt="zonda"/>

## General Information

Zonda is an web interface that facilitates generating EXTPAR data on ICON triangular grids for
reasearch and on-demand simulations.

## Main Features

- Users can configure the desired grid and EXTPAR options in the web UI.
- Users will get a preview of the domain on a 3D globe or 2D map.
- Up to four (nested) domains can be generated at once.
- ICON grid and external parameters are provided to the user via download link.
- Data processing takes place at ETH Zurich servers.
- A GitHub account is required to submit requests.

## Web Interface

The graphical web interface consists of the following elements:

- Title and button to switch between light and dark mode
- Domain Specifications
- EXTPAR Parameters
- JSON Request for GitHub

## Quick Start

!!! tip
    For all variables that can be set by the user, there is a tooltip available by
    hovering the :material-information-outline:{ title="I am a tooltip" } icon.


### 1. Defining the base grid

Under **RnBk Specifications**, define the grid values `n` (`grid_root`) and `k` (`grid_level`)
for the base grid and adapt the output file name.

Below, estimations of horizontal resolution and number of grid points are displayed for the
base domain, based on the `RnBk` values.

### 2. Setting up domains

**Domain 1** displays all settings for the base domain, which can be a local or a global one
(see the toggle below). 

You can add further domains by clicking on the :octicons-plus-circle-16: icon on the right.
Note that the variable `parent_id` refers to the domain number, e.g., if Domain 2 should be
a nest within Domain 1, set `parent_id` to 1 for Domain 2.

You can delete domains by clicking on the :fontawesome-regular-trash-can: icon.

### 3. Preview

Once valid grid numbers are put in, your domain(s) is/are shown in the preview on the left-hand-side.
You can toggle between the 3D (globe) and 2D (map) views.

### 4. EXTPAR Parameters

These are all the relevant variables required by EXTPAR. For more information on these variables,
click on the variable name that leads to the relevant section of the documentation.

### 5. JSON Request

The web interface automatically generates a code snippet in the popular
[JSON format :material-open-in-new:](https://en.wikipedia.org/wiki/JSON){:target="_blank"} based on the values entered by the user in
the sections above.

If everything is alright, use the "Click to copy" button and then "Open Issue", which redirects
you to the data request issue on GitHub. There, you need to paste the JSON code (which has been
automatically copied to the clipboard) by replacing the line `PASTE_YOUR_REQUEST_HERE`.
For more information about the data request and to check if your JSON is correctly inserted,
click on "Preview" in the issue box before you submit it.

!!! note
    The [zonda-request :material-open-in-new:](https://github.com/C2SM/zonda-request){:target="_blank"} repository is publicly available. 
    However, to create a new 
    [data request issue :material-open-in-new:](https://github.com/C2SM/zonda-request/issues/new?template=data-request.md){:target="_blank"}
    (and therefore submit a request), you need to have a [GitHub user account :material-open-in-new:](https://github.com/signup){:target="_blank"}.