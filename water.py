import requests
import json
import os
import geopandas as gpd
import requests
import matplotlib.pyplot as plt
import contextily as ctx
from shapely.geometry import Point
import numpy as np
from matplotlib.colors import ListedColormap
import pyproj
from shapely.geometry import Polygon, LineString, MultiPolygon, shape
from shapely.geometry import box
from matplotlib.widgets import Button
from shapely.ops import unary_union
from scipy.ndimage import geometric_transform
from skimage.transform import warp_polar
import cv2

def fetch_water_bodies(south, west, north, east, filename="overpass_data.json"):
    """
    Fetch water bodies within the bounding box using Overpass API.
    If the data is already stored locally, it loads from the file.
    Otherwise, it fetches from the Overpass API and stores the result in a file.
    Parameters:
        - south, west, north, east (float): Bounding box coordinates.
        - filename (str): Name of the file to store/fetch local data.
    Returns:
        - dict: JSON response from the Overpass API or local file.
    """
    
    # If file exists, load data from it
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            return json.load(f)
    
    # Define a bounding box string
    bbox = f"{south},{west},{north},{east}"
    
    # Overpass QL for fetching water bodies within or intersecting the bounding box
    overpass_ql = f"""
    [out:json][timeout:25];
    (
        way["water"="lake"]({bbox});
        way["natural"="coastline"]({bbox});
        way["waterway"="riverbank"]({bbox});
        // query part for: “natural=water and type=multipolygon”
        node["natural"="water"]["type"="multipolygon"]({bbox});
        way["natural"="water"]["type"="multipolygon"]({bbox});
        relation["natural"="water"]["type"="multipolygon"]({bbox});
        // query part for: “natural=water”
        node["natural"="water"]({bbox});
        way["natural"="water"]({bbox});
        relation["natural"="water"]({bbox});
    );
    out body;
    >;
    out skel qt;
    """
    
    url = "https://overpass-api.de/api/interpreter"
    response = requests.get(url, params={"data": overpass_ql})
    response.raise_for_status()
    
    # Save the data to a file
    with open(filename, 'w') as f:
        json.dump(response.json(), f)
    
    return response.json()

def overpass_json_to_geodf(data):
    # Convert nodes and ways from the JSON into a list of shapely geometries
    nodes = {node["id"]: (node["lon"], node["lat"]) for node in data["elements"] if node["type"] == "node"}
    geometries = []
    for element in data["elements"]:
        if element["type"] == "way":
            try:
                way_nodes = [nodes[node] for node in element["nodes"]]
                # If the way is not closed, add the first node at the end to close it
                if way_nodes[0] != way_nodes[-1]:
                    way_nodes.append(way_nodes[0])
                polygon = shape({
                    "type": "Polygon",
                    "coordinates": [way_nodes]
                })
                # Check if the polygon is valid
                if not polygon.is_valid:
                    # Attempt to fix invalid polygons
                    polygon = polygon.buffer(0)
                geometries.append(polygon)
            except Exception as e:
                print(str(e))
                continue  # Just skip any problematic ways

    # Use unary_union to merge overlapping or adjacent polygons
    merged_geometries = unary_union(geometries)

    # Create a GeoDataFrame from the merged geometries
    # If unary_union results in a single Polygon, we wrap it in a list to create a GeoDataFrame
    if isinstance(merged_geometries, (Polygon, MultiPolygon)):
        gdf = gpd.GeoDataFrame(geometry=[merged_geometries], crs="EPSG:4326")
    else:
        # Handle other types (e.g., if there's a mix of geometry types)
        gdf = gpd.GeoDataFrame(geometry=list(merged_geometries), crs="EPSG:4326")

    return gdf


vertices = []
globalWest = 0
globalSouth = 0
globalGrid = None
globalResolution = 0
globalExtentMap = []

def plot_with_osm_background_and_grid(gdf, grid, south, west, north, east, resolution=0.00005):
    global globalGrid
    global globalWest
    global globalSouth
    global globalResolution
    globalResolution = resolution
    globalGrid = grid

    # Set the initial CRS to WGS 84 if it's not set already
    if gdf.crs is None:
        gdf = gdf.set_crs(epsg=4326)

    # Ensure the GeoDataFrame uses the correct coordinate system
    gdf = gdf.to_crs(epsg=3857)

    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(12, 12))

    # Plot the polygons on the axis
    #gdf.plot(ax=ax, alpha=0.5, color='blue', edgecolor='k', zorder=2)

    # Adjust the grid extent to the map's coordinate system (from lat-lon to Web Mercator)
    transformer = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)

    new_west, new_south = transformer.transform(west, south)
    globalWest = west
    globalSouth = south
    new_east, new_north = transformer.transform(east, north)

    extent_grid = [new_west, new_east, new_south, new_north]
    global globalExtentMap
    globalExtentMap = extent_grid

    # Create a custom colormap where '0' values (unknown) are set as transparent
    colors = [(0, 0, 0, 0), (0, 0, 0, 1)]  # RGBA colors (transparent, black)
    cmap = ListedColormap(colors, N=2)

    # Draw the grid as an image with the custom colormap
    ax.imshow(grid, cmap=cmap, extent=extent_grid, origin='lower', zorder=3, aspect='auto')

    # Add OSM tiles as background
    ctx.add_basemap(ax = ax, source=ctx.providers.OpenStreetMap.Mapnik, zorder=1)

    def onclick(event):
        print("Got click!")
        toolbar = plt.get_current_fig_manager().toolbar
        if not (event.inaxes == ax and toolbar.mode == ''):
            return
        global vertices
        if event.inaxes != ax:
            return

        # Get the click coordinates in Web Mercator
        click_x, click_y = event.xdata, event.ydata

        transformer_reverse = pyproj.Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)
        # Transform these coordinates back to lat/lon
        lon, lat = transformer_reverse.transform(click_x, click_y)

        global globalResolution
        ix = int((lon - globalWest) / globalResolution)
        iy = int((lat - globalSouth) / globalResolution)
        vertices.append((ix, iy))
        
        # Plot the point for visual confirmation
        plt.plot(click_x, click_y, 'ro')
        plt.draw()

    def onkeypress(event):
        """Handle key press events."""
        global globalGrid, vertices
        if event.key == 't':
            if vertices:
                # Transform vertices from plot pixel coordinates to grid indices
                transformed_vertices = [(int(x), int(globalGrid.shape[0] - 1 - y)) for x, y in vertices]
                new_water_area = Polygon(transformed_vertices)
                for ix in range(globalGrid.shape[1]):
                    for iy in range(globalGrid.shape[0]):
                        # Check contains using grid indices
                        if new_water_area.contains(Point(ix, grid.shape[0] - 1 - iy)):
                            globalGrid[iy, ix] = 1
                np.save("updated_grid.npy", globalGrid)
                grid_binary = (grid * 255).astype(np.uint8)
                cv2.imwrite("updated_image.png", grid_binary)
                vertices = []  # Clear vertices for next polygon
                ax.clear()
                colors = [(0, 0, 0, 0), (0, 0, 0, 1)]  # RGBA colors (transparent, black)
                cmap = ListedColormap(colors, N=2)
                global globalExtentMap
                ax.imshow(globalGrid, cmap=cmap, extent=globalExtentMap, origin='lower', zorder=3, aspect='auto')
                ctx.add_basemap(ax = ax, source=ctx.providers.OpenStreetMap.Mapnik, zorder=1)
                
                plt.draw()

        cid = fig.canvas.mpl_connect('button_press_event', onclick)

    fig.canvas.mpl_connect('key_press_event', onkeypress)
    fig.canvas.mpl_connect('button_press_event', onclick)

    # Show the plot
    plt.show()

def create_occupancy_grid(data, south, west, north, east, resolution=0.00005, filename = "grid_binary.npy"):
    """
    Create an occupancy grid from water body data, checking for cell intersection with geometries.
    Parameters:
        - data: GeoDataFrame containing the water bodies.
        - south, west, north, east (float): Bounding box coordinates.
        - resolution (float): Size of each cell in lat/lon degrees.
    Returns:
        - numpy.ndarray: Occupancy grid.
    """ 
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            grid = np.load(filename)
            grid_binary = (grid * 255).astype(np.uint8)
            cv2.imwrite("image.png", grid_binary)
            return grid
    # Calculate the number of cells in the grid
    grid_x = int((east - west) / resolution) + 1
    grid_y = int((north - south) / resolution) + 1
    grid = np.zeros((grid_y, grid_x))
    
    # Iterate through each polygon in the data
    # for geometry in data['geometry']:
    #     if geometry is not None and not geometry.is_empty:
    #         for ix in range(grid_x):
    #             for iy in range(grid_y):
    #                 # Calculate the bounds of the current cell
    #                 x_min = west + ix * resolution
    #                 y_min = south + iy * resolution
    #                 x_max = x_min + resolution
    #                 y_max = y_min + resolution

    #                 # Create a box for the current cell
    #                 cell_box = box(x_min, y_min, x_max, y_max)

    #                 # Check for intersection between the cell and the geometry
    #                 if geometry.intersects(cell_box):
    #                     grid[iy, ix] = 1  # Mark as occupied
                        
    # Optional: Mark free cells (cells not occupied by water)
    # grid[grid == 0] = -1
    
    # Save the data to a file
    np.save(filename, grid)
    grid_binary = (grid * 255).astype(np.uint8)
    cv2.imwrite("image.png", grid_binary)
    return grid

def fix_geometry(geometry):
    """
    Attempt to fix an open geometry by ensuring it is properly closed.
    This function handles LineString and Polygon geometries.
    """
    if isinstance(geometry, LineString):
        # Convert a LineString to a Polygon if it's not closed
        if not geometry.is_ring:
            # Attempt to close the LineString by creating a Polygon
            return Polygon(geometry)
    elif isinstance(geometry, Polygon):
        # Ensure the Polygon is properly closed
        exterior = geometry.exterior
        if not exterior.is_ring:
            # Close the polygon by creating a new one with the closed exterior
            return Polygon(list(exterior.coords) + [exterior.coords[0]])
    return geometry

def main():
    # Define a bounding box (e.g., around a part of a city or lake)
    #south, west, north, east = 42.838070, -71.006272, 42.862226, -70.969131  # Bounding box around lake attitash
    #south, west, north, east = 42.274460, -71.759992, 42.280359, -71.752954 # Bounding box around lake Quinsigamond
    #south, west, north, east = 42.272388, -71.807720, 42.274761, -71.803917 # Bounding box around the lab
    south, west, north, east = 42.02697586726279, -71.86015011725914, 42.04778634200546, -71.830713667753 # Webster lake
    data = fetch_water_bodies(south, west, north, east)
    
    # Convert Overpass JSON data to GeoDataFrame
    water_gdf = overpass_json_to_geodf(data)

    if water_gdf.empty:
        print("No water polygons found in the specified bounding box.")
    else:
        grid = create_occupancy_grid(water_gdf, south, west, north, east)
        plot_with_osm_background_and_grid(water_gdf, grid, south, west, north, east)

if __name__ == "__main__":
    main()
