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
                # Checking if the way is closed to determine if it's a polygon
                if way_nodes[0] == way_nodes[-1]:
                    geometries.append({
                        "type": "Polygon",
                        "coordinates": [way_nodes]
                    })
                else:
                    geometries.append({
                        "type": "LineString",
                        "coordinates": way_nodes
                    })
            except:
                continue  # Just skip any problematic ways

    # Convert geometries list into a GeoDataFrame
    gdf = gpd.GeoDataFrame.from_features({
        "type": "FeatureCollection", 
        "features": [{"type": "Feature", "geometry": geom, "properties": {}} for geom in geometries]  # Added "properties": {}
    })
    return gdf

def plot_with_osm_background_and_grid(gdf, grid, south, west, north, east, resolution=0.0001):
    # Set the initial CRS to WGS 84 if it's not set already
    if gdf.crs is None:
        gdf = gdf.set_crs(epsg=4326)

    # Ensure the GeoDataFrame uses the correct coordinate system
    gdf = gdf.to_crs(epsg=3857)

    # Create a figure and axis object
    fig, ax = plt.subplots(figsize=(12, 12))

    # Plot the polygons on the axis
    gdf.plot(ax=ax, alpha=0.5, color='blue', edgecolor='k', zorder=2)

    # Determine the extent in Web Mercator coordinates
    x_min, y_min, x_max, y_max = ax.axis()
    extent_map = [x_min, x_max, y_min, y_max]

    # Adjust the grid extent to the map's coordinate system (from lat-lon to Web Mercator)
    transformer = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)

    new_west, new_south = transformer.transform(west, south)
    new_east, new_north = transformer.transform(east, north)

    extent_grid = [new_west, new_east, new_south, new_north]

    # Create a custom colormap where '0' values (unknown) are set as transparent
    colors = [(0, 0, 0, 0), (0, 0, 0, 1)]  # RGBA colors (transparent, black)
    cmap = ListedColormap(colors, N=2)

    # Draw the grid as an image with the custom colormap
    ax.imshow(grid, cmap=cmap, extent=extent_grid, origin='lower', zorder=3, aspect='auto')

    # Add OSM tiles as background
    ctx.add_basemap(ax = ax, source=ctx.providers.OpenStreetMap.Mapnik, zorder=1)

    # Show the plot
    plt.show()

def create_occupancy_grid(data, south, west, north, east, resolution=0.0001, filename = "grid_binary.npy"):
    """
    Create an occupancy grid from water body data.
    Parameters:
        - data: GeoDataFrame containing the water bodies.
        - south, west, north, east (float): Bounding box coordinates.
        - resolution (float): Size of each cell in lat/lon degrees.
    Returns:
        - numpy.ndarray: Occupancy grid.
    """
    
    if os.path.exists(filename):
        with open(filename, 'r') as f:
            return np.load(filename)

    # 0 for unknown, 1 for occupied, -1 for free
    grid_x = int((east - west) / resolution)+1
    grid_y = int((north - south) / resolution)+1
    grid = np.zeros((grid_y, grid_x))
    
    # Iterate through each polygon in the data
    for geometry in data['geometry']:
        if geometry is not None and not geometry.is_empty:
            for x in np.arange(west, east, resolution):
                for y in np.arange(south, north, resolution):
                    cell_center = Point(x + resolution/2, y + resolution/2)
                    if geometry.contains(cell_center):
                        ix = int((x - west) / resolution)
                        iy = int((y - south) / resolution)
                        grid[iy, ix] = 1  # Mark as occupied
                        
    # Optional: Mark free cells (cells not occupied by water)
    # grid[grid == 0] = -1
    # Save the data to a file
    np.save(filename, grid)
    return grid

def main():
    # Define a bounding box (e.g., around a part of a city or lake)
    south, west, north, east = 42.838070, -71.006272, 42.862226, -70.969131  # Bounding box around lake attitash
    data = fetch_water_bodies(south, west, north, east)

    # Now, 'data' contains the water body polygons in that bounding box.
    # You can further process this data or visualize it using a tool or library

    # Convert Overpass JSON data to GeoDataFrame
    water_gdf = overpass_json_to_geodf(data)

    if water_gdf.empty:
        print("No water polygons found in the specified bounding box.")
    else:
        grid = create_occupancy_grid(water_gdf, south, west, north, east)
        plt.imshow(grid, cmap='gray', origin='lower')
        plot_with_osm_background_and_grid(water_gdf, grid, south, west, north, east)

if __name__ == "__main__":
    main()
