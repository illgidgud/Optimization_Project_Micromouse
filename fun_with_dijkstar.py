import sys
# Replace this with whichever path your package is, if not already included
# in sys.path
sys.path.append("/Library/Frameworks/Python.framework/Versions/3.10/lib/python3.10/site-packages")
from pathlib import Path
from dijkstar import Graph, find_path
import json
import time
import math
from functools import cmp_to_key
import numpy as np
from matplotlib import pyplot as plt

# Slopes generator
# Reduce a slope
def reduce_tuple(iter):
    """
    Reduce the iterable (a, b)
    """
    a, b = iter[0], iter[1]
    d = math.gcd(a, b)
    return (int(a/d), int(b/d))

# Get available slopes from a specified position
# Note: This method uses a set, which means order will change each run
def get_slopes(current_col, current_row, row, col):
    """
    Return list of slopes for point ('i', 'j') in a
    'row' * 'col' grid. A slope is represented as an
    irreducible list [a, b], in the right-upward direction (i.e b > 0,
    or b = 0 and a > 0)
    """
    slope_set = set()
    if current_col < col:
        slope_set.add((1, 0))
    for row_index in range(current_row + 1, row + 1):
        for col_index in range(1, col + 1):
            slope_set.add(reduce_tuple((col_index - current_col, row_index - current_row)))
    return list(slope_set)

# Generate slopes for square grids, put save = "yes" for local save of slopes.
# Return None if save = "yes"
def generate_slopes(size, save="no"):
    """
    Generate slopes for square grids, put save = "yes" for local save of slopes.
    Return None if save = "yes"
    """
    slope_dict = dict()
    for current_row in range (1, size + 1):
        for current_col in range (1, size + 1):
            slope_dict[f"{current_col}_{current_row}"] = get_slopes(
                current_col, current_row, size, size,
            )
    if save == "yes":
        path = Path(__file__).parent/"Slopes"/f"Size{size}"
        with open(path, "w") as f:
            json.dump(slope_dict, f, indent=1)
            return None
    return slope_dict

# Split nodes_info keys, for both models
def split_dict_key(key):
    """
    Convert dict key 'i_j' to list [i, j]
    """
    coords = key.split("_")
    return [int(coords[0]), int(coords[1])]

class Constant:
    """
    Turning time: After normalizing grid size to 1 and velcity to 1.
    Additional assumption: Turning time for acute angles 
    is twice as long as that for obtuse angles.
    """
    # Based on information in the video "The Fastest Maze-Solving Competition On Earth"
    # of channel Veritasium, from 10:40 to 11:20.
    acute =  7.5
    right_or_obtuse =  3.75

# For model 2
def get_angle_cost_for_model_2(slope1, slope2):
    """
    return the cost between the two slopes (in their positive direction)    
    """
    if slope1[0] * slope2[1] - slope1[1] * slope2[0] == 0: return 0
    else:
        inner_product = slope1[0] * slope2[0] + slope1[1] * slope2[1]
        if inner_product < 0: return Constant.acute
        else: return Constant.right_or_obtuse

def cost_function(u, v, edge, prev_edge):
    current_length, current_slope = edge
    if prev_edge:
        prev_slope = prev_edge[1]
    else:
        return current_length
    return current_length + get_angle_cost_for_model_2(prev_slope, current_slope)

"""
Tổng kết: quá trình thực hiện cho model 1
1. Đọc thông tin của lưới ô vuông
2. (Tuỳ chọn) Xoá đi các đỉnh là đầu mút của tường. Chú ý thiết kế code
sao cho việc xoá bằng list.remove() không trả về lỗi khi không có phần
tử cần xoá.
3. Tạo ra dữ liệu về nodes_info (dictionary) (gồm các cặp đỉnh và hệ số góc tối giản
có thể đi).
4. Từ dữ liệu nodes_info trên, tiếp tục xoá đi các đỉnh "không quan trọng",
là các đỉnh mà chẳng có một hướng nào đi được.
5. Từ active_nodes_info, tạo ra các đỉnh của đồ thị và một quan hệ thứ tự
để phiên dịch từ các đỉnh trừu tượng ra con số và ngược lại.
6. Thêm các cạnh. Có hai kiểu cạnh:
- Một là cạnh giữa hai "nhánh" khác nhau,
- Hai là cạnh giữa hai đỉnh kề nhau trong cùng một nhánh.
7. Giải đồ thị, in ra kết quả về các đỉnh. Lưu ý, có thể có sự trùng lặp
về đỉnh trong lời giải này (nằm ở các nhánh khác nhau, ứng với các hướng
quay khác nhau). Cần lọc các đỉnh lặp lại để thu được dãy đỉnh là kết quả
cuối cùng. Có thể tích hợp luôn phần vẽ kết quả cuối cùng, hoặc lưu ra một file rồi
vẽ sau.
"""

# Get Euclidean distance, with chosen round-digit
def get_Cartesian_length(slope, round_digit=5):
    """
    return Cartesian length of slope
    """
    a, b = slope[0], slope[1]
    return round(math.sqrt(a*a + b*b), round_digit)

# Get angle cost
def get_angle_cost(slope1, slope2):
    """
    return the cost between the two slopes (in their positive direction)    
    """
    inner_product = slope1[0] * slope2[0] + slope1[1] * slope2[1]
    if inner_product > 0: return [Constant.acute, Constant.right_or_obtuse]
    elif inner_product == 0: return [Constant.right_or_obtuse, Constant.right_or_obtuse]
    else: return [Constant.right_or_obtuse, Constant.acute]

# Check intersection
def check_intersection(point1, point2, edge):
    """
    Check intersection of segment ['point1', 'point2'] with 'edge'
    """
    a1, b1, a2, b2 = edge[0][0], edge[0][1], edge[1][0], edge[1][1]
    x1, y1, x2, y2 = point1[0], point1[1], point2[0], point2[1]
    """
    Solve the following system
    u*(a1 - a2) - v*(x1 - x2) = x2 - a2
    u*(b1 - b2) - v*(y1 - y2) = y2 - b2
    Condition of intersection: u, v in [0, 1]
    """
    d = (a2 - a1)*(y1 - y2) + (b1 - b2)*(x1 - x2)
    d1 = (a2 - x2)*(y1 - y2) + (y2 - b2)*(x1 - x2)
    d2 = (a1 - a2)*(y2 - b2) + (b2 - b1)*(x2 - a2)
    if d == 0:
        if d1 != 0 or d2 != 0: return False
        else:
            a12 = a1 - a2
            if a12 == 0:
                b12 = b1 - b2
                y12 = y1 - y2
                y1b2 = y1 - b2
                y2b1 = y2 - b1
                y1b1 = y1 - b1
                y2b2 = y2 - b2
                if b12 > 0 and y12 > 0:
                    return (y1b2 >= 0 and y2b1 <= 0)
                elif b12 > 0 and y12 < 0:
                    return (y2b2 >= 0 and y1b1 <= 0)
                elif b12 < 0 and y12 > 0:
                    return (y2b2 <= 0 and y1b1 >= 0)
                else:
                    return (y1b2 <= 0 and y2b1 >= 0)
            else:
                x12 = x1 - x2
                x1a2 = x1 - a2
                x2a1 = x2 - a1
                x1a1 = x1 - a1
                x2a2 = x2 - a2
                if a12 > 0 and x12 > 0:
                    return (x1a2 >= 0 and x2a1 <= 0)
                elif a12 > 0 and x12 < 0:
                    return (x2a2 >= 0 and x1a1 <= 0)
                elif a12 < 0 and x12 > 0:
                    return (x2a2 <= 0 and x1a1 >= 0)
                else:
                    return (x1a2 <= 0 and x2a1 >= 0)
                
    else:
        cond_u = (all([d > 0, d1 >= 0, d1 - d <= 0]) or all([d < 0, d1 <= 0, d1 - d >= 0]))
        cond_v = (all([d > 0, d2 >= 0, d2 - d <= 0]) or all([d < 0, d2 <= 0, d2 - d >= 0]))
        return all([cond_u, cond_v])

# Get reach from a point
def get_furthest_reach(point, slope, row, col, edge_list):
    """
    Get furthest expansion from 'point' in the direction given by 'slope',
    in a 'row'*'col' grid with 'edge_list' as obstacles.
    """
    current_col, current_row = point[0], point[1]
    col_slope, row_slope = slope[0], slope[1]
    # Get maximum column reach
    if col_slope == 0: max_col_reach = col
    elif col_slope > 0: max_col_reach = int((col - current_col) / col_slope) + 1
    else: max_col_reach = int((current_col - 1) / (-col_slope)) + 1
    # Get maximum row reach
    if row_slope == 0: max_row_reach = row
    else: max_row_reach = int((row - current_row) / row_slope) + 1

    min_reach = 0
    max_reach = min(max_row_reach, max_col_reach)
    while max_reach - min_reach > 1:
        temp = min_reach
        new_reach = int((max_reach + min_reach) / 2)
        min_reach = new_reach
        test_point = [current_col + new_reach * col_slope, current_row + new_reach * row_slope]
        for edge in edge_list:
            if check_intersection(point, test_point, edge):
                max_reach = new_reach
                min_reach = temp
                break
    return min_reach

# Sort the edge_list, to speed up categorization
def sort_edge_list(edge_list):
    """
    Sort edge_list by length of edge, from longest to shortest, using
    taxicab distance
    """
    return sorted(
        edge_list, 
        key=lambda x: abs(x[0][0] - x[1][0]) + abs(x[0][1] - x[1][1]),
        reverse=True,
    )

# Convert grid coordinates to num, and vice versa. Used in nodes dictionary
# for graph generation and translation.
def coords_to_num(current_col, current_row, col):
    """
    Note: no row index needed
    """
    return col*(current_row - 1) + current_col

def num_to_coords(num, col):
    """
    Note: No row index needed
    """
    if num / col == int(num / col):
        current_row = num / col
    else: 
        current_row = int(num / col) + 1
    current_col = num - col*(current_row - 1)
    return [current_col, current_row]

# Get points on a segment, for first-stage filtering of redundant vertices,
# among other purposes.
def get_points_between(point1, point2):
    """
    Return list of all integral points between 'point1' and 'point2'
    """
    a1, b1, a2, b2 = point1[0], point1[1], point2[0], point2[1]
    a, b = a2 - a1, b2 - b1
    if b > 0 or (b == 0 and a > 0):
        d = math.gcd(a, b)
        i, j = int(a/d), int(b/d)
        return [[a1 + num * i, b1 + num * j] for num in range(0, d + 1)]
    else:
        a, b = -a, -b
        d = math.gcd(a, b)
        i, j = int(a/d), int(b/d)
        return [[a2 + num * i, b2 + num * j] for num in range(0, d + 1)]

# Also for filtering purpose
def is_empty_dict(dict):
    """
    Checks for empty dictionaries (values are empty lists)
    """
    return all(len(dict[key]) == 0 for key in dict.keys())

# Comparing nodes, to generate nodes dictionary
def compare(node1, node2):
    """
    Compare nodes of the following structure
    node = "{col}_{row}_{slope_col}_{slope_row}_direction"
    return -1 if node1 < node2
    """
    list1 = [int(char) for char in node1.split("_")]
    list2 = [int(char) for char in node2.split("_")]
    if list1 < list2: return -1
    elif list1 == list2: return 0
    else: return 1

def solve_with_first_model(size, index, visualize = "no"):
    start_time = time.time()
    # List of all points (to be filtered and used later)
    # Note: j before i
    points_list = [[i, j] for j in range(1, size + 1) for i in range(1, size + 1)]

    nodes_info = {
        f"{i}_{j}": {
            "begin": [], 
            "middle": [], 
            "end": [],
            "reach": [],
        } 
        for j in range(1, size + 1) for i in range(1, size + 1)
    }

    # Two options for nodes_slopes, one through local save, another through
    # in-program generation. Testing seems to prefer the save and read option.
    slopes_path = Path(__file__).parent/"Slopes"/f"Size{size}"
    with open(slopes_path, "r") as f:
        node_slopes = json.load(f)

    # node_slopes = generate_slopes(size=size)

    grid_path = Path(__file__).parent/"Samples"/f"Size{size}"/f"sample{index}.json"
    with open(grid_path, "r") as f:
        grid_info = json.load(f)

    # Get edges, then sort (for faster execution)
    edges = grid_info["edges"]
    edges = sort_edge_list(edge_list=edges)

    # Filter redundant points from points_list
    for edge in edges:
        points_to_remove = get_points_between(edge[0], edge[1])
        for point in points_to_remove:
            # remove from points_list
            try:
                points_list.remove(point)
            except:
                pass

            # remove from nodes_info
            key = f"{point[0]}_{point[1]}"
            try:
                del nodes_info[key]
            except:
                pass
    
    # First validity check
    start = grid_info["start"]
    target = grid_info["target"]
    if start not in points_list: 
        print(f"Invalid: Start is on the wall")
        return None
    if target not in points_list: 
        print(f"Invalid: Target is on the wall")
        return None

    # Nodes_info
    for current_col, current_row in points_list:
        point = [current_col, current_row]
        for slope in node_slopes[f"{current_col}_{current_row}"]:
            max_reach = get_furthest_reach(point, slope, size, size, edges)
            if max_reach > 0:
                nodes_info[f"{current_col}_{current_row}"]["begin"].append(slope)
                nodes_info[f"{current_col}_{current_row}"]["reach"].append(max_reach)
                for num in range(1, max_reach):
                    new_col = current_col + num * slope[0]
                    new_row = current_row + num * slope[1]
                    try:
                        node_slopes[f"{new_col}_{new_row}"].remove(slope)
                    except: 
                        pass
                    nodes_info[f"{new_col}_{new_row}"]["middle"].append(slope)
                new_col = current_col + max_reach * slope[0]
                new_row = current_row + max_reach * slope[1]
                try:
                    node_slopes[f"{new_col}_{new_row}"].remove(slope)
                except:
                    pass
                nodes_info[f"{new_col}_{new_row}"]["end"].append(slope)

    # Second validity check (start and target both in nodes_info)
    if is_empty_dict(nodes_info[f"{start[0]}_{start[1]}"]):
        print(f"Invalid: Start is disconnected")
        return None

    if is_empty_dict(nodes_info[f"{target[0]}_{target[1]}"]):
        print(f"Invalid: Target is disconnected")
        return None

    else:
        # Create graph
        graph = Graph()

        # Further filter unreachable nodes from list
        for current_col, current_row in points_list:
            key = f"{current_col}_{current_row}"
            if is_empty_dict(nodes_info[key]):
                del nodes_info[key]

        # Create nums_to_nodes and nodes_to_nums dictionaries, for graphing
        abstract_nodes_list = []
        for key, slopes in nodes_info.items():
            for slope in slopes["begin"]:
                abstract_nodes_list.append(key + f"_{slope[0]}_{slope[1]}_1")
                abstract_nodes_list.append(key + f"_{slope[0]}_{slope[1]}_-1")
            for slope in slopes["middle"]:
                abstract_nodes_list.append(key + f"_{slope[0]}_{slope[1]}_1")
                abstract_nodes_list.append(key + f"_{slope[0]}_{slope[1]}_-1")
            for slope in slopes["end"]:
                abstract_nodes_list.append(key + f"_{slope[0]}_{slope[1]}_1")
                abstract_nodes_list.append(key + f"_{slope[0]}_{slope[1]}_-1")

        abstract_nodes_list = sorted(abstract_nodes_list, key=cmp_to_key(compare))
        nodes_to_nums = dict()
        nums_to_nodes = dict()
        for i in range(len(abstract_nodes_list)):
            nodes_to_nums[abstract_nodes_list[i]] = i
            nums_to_nodes[i] = abstract_nodes_list[i]

        # Create edges
        # Crate adjacent edges first, note that I pop the start and target
        # points later from nodes_info
        for key in nodes_info.keys():
            point = split_dict_key(key)
            length = len(nodes_info[key]["begin"])
            for i in range(length):
                slope = nodes_info[key]["begin"][i]
                reach = nodes_info[key]["reach"][i]
                slope_key = f"_{slope[0]}_{slope[1]}_"
                unit_length = get_Cartesian_length(slope)
                for j in range(reach):
                    first_point_key = f"{point[0] + j * slope[0]}_{point[1] + j * slope[1]}" + slope_key
                    first_node_pos_key = nodes_to_nums[first_point_key + "1"]
                    first_node_neg_key = nodes_to_nums[first_point_key + "-1"]

                    second_point_key = f"{point[0] + (j + 1) * slope[0]}_{point[1] + (j + 1) * slope[1]}" + slope_key
                    second_node_pos_key = nodes_to_nums[second_point_key + "1"]
                    second_node_neg_key = nodes_to_nums[second_point_key + "-1"]

                    # Add edges
                    graph.add_edge(first_node_pos_key, second_node_pos_key, unit_length)
                    graph.add_edge(second_node_neg_key, first_node_neg_key, unit_length)

        # Add edges between the same nodes with different angles
        # Get start and target
        start_key = f"{start[0]}_{start[1]}"
        target_key = f"{target[0]}_{target[1]}"
        start_slopes = nodes_info.pop(start_key)
        target_slopes = nodes_info.pop(target_key)

        start_slope_list = start_slopes["begin"] + start_slopes["middle"] + start_slopes["end"]
        target_slope_list = target_slopes["begin"] + target_slopes["middle"] + target_slopes["end"]

        # Add start and target to edge_dict
        first_start = start_slope_list.pop(0)
        first_start_key = start_key + f"_{first_start[0]}_{first_start[1]}_"
        first_start_pos_key = nodes_to_nums[first_start_key + "1"]
        first_start_neg_key = nodes_to_nums[first_start_key + "-1"]

        first_target = target_slope_list.pop(0)
        first_target_key = target_key + f"_{first_target[0]}_{first_target[1]}_"
        first_target_pos_key = nodes_to_nums[first_target_key + "1"]
        first_target_neg_key = nodes_to_nums[first_target_key + "-1"]

        # start and target: all edges are 0
        for second_start in start_slope_list:
            second_start_key = start_key + f"_{second_start[0]}_{second_start[1]}_"
            second_start_pos_key = nodes_to_nums[second_start_key + "1"]
            second_start_neg_key = nodes_to_nums[second_start_key + "-1"]
            # Add edges
            graph.add_edge(first_start_pos_key, second_start_neg_key, 0)
            graph.add_edge(first_start_pos_key, second_start_neg_key, 0)
            graph.add_edge(second_start_neg_key, first_start_pos_key, 0)
            graph.add_edge(first_start_neg_key, second_start_pos_key, 0)
            graph.add_edge(second_start_pos_key, first_start_neg_key, 0)
            graph.add_edge(first_start_pos_key, second_start_pos_key, 0)
            graph.add_edge(second_start_pos_key, first_start_pos_key, 0)
            graph.add_edge(first_start_neg_key, second_start_neg_key, 0)
            graph.add_edge(second_start_neg_key, first_start_neg_key, 0)

        for second_target in target_slope_list:
            second_target_key = target_key + f"_{second_target[0]}_{second_target[1]}_"
            second_target_pos_key = nodes_to_nums[second_target_key + "1"]
            second_target_neg_key = nodes_to_nums[second_target_key + "-1"]
            # Add edges
            graph.add_edge(first_target_pos_key, second_target_neg_key, 0)
            graph.add_edge(second_target_neg_key, first_target_pos_key, 0)
            graph.add_edge(first_target_neg_key, second_target_pos_key, 0)
            graph.add_edge(second_target_pos_key, first_target_neg_key, 0)
            graph.add_edge(first_target_pos_key, second_target_pos_key, 0)
            graph.add_edge(second_target_pos_key, first_target_pos_key, 0)
            graph.add_edge(first_target_neg_key, second_target_neg_key, 0)
            graph.add_edge(second_target_neg_key, first_target_neg_key, 0)

        for key, slopes in nodes_info.items():
            begin_slopes = slopes["begin"]
            middle_slopes = slopes["middle"]
            end_slopes = slopes["end"]
            begin_slope_length = len(begin_slopes)
            middle_slopes_length = len(middle_slopes)
            end_slopes_length = len(end_slopes)

            # First loop
            for i in range(begin_slope_length - 1):
                first_slope = slopes["begin"][i]
                first_slope_key = key + f"_{first_slope[0]}_{first_slope[1]}_"
                first_node_pos_key = nodes_to_nums[first_slope_key + "1"]
                first_node_neg_key = nodes_to_nums[first_slope_key + "-1"]
                for j in range(i + 1, begin_slope_length):
                    #Initiate materials
                    second_slope = slopes["begin"][j]
                    angle_cost = get_angle_cost(first_slope, second_slope)[0]
                    second_slope_key = key + f"_{second_slope[0]}_{second_slope[1]}_"
                    second_node_pos_key = nodes_to_nums[second_slope_key + "1"]
                    second_node_neg_key = nodes_to_nums[second_slope_key + "-1"]

                    # Add edges
                    graph.add_edge(second_node_neg_key, first_node_pos_key, angle_cost)
                    graph.add_edge(first_node_neg_key, second_node_pos_key, angle_cost)

                for second_slope in middle_slopes:
                    costs = get_angle_cost(first_slope, second_slope)
                    angle_cost = costs[0]
                    other_angle_cost = costs[1]
                    second_slope_key = key + f"_{second_slope[0]}_{second_slope[1]}_"
                    second_node_pos_key = nodes_to_nums[second_slope_key + "1"]
                    second_node_neg_key = nodes_to_nums[second_slope_key + "-1"]

                    # Add edges
                    graph.add_edge(second_node_neg_key, first_node_pos_key, angle_cost)
                    graph.add_edge(first_node_neg_key, second_node_pos_key, angle_cost)

                    graph.add_edge(second_node_pos_key, first_node_pos_key, other_angle_cost)
                    graph.add_edge(first_node_neg_key, second_node_neg_key, other_angle_cost)

                for second_slope in end_slopes:
                    other_angle_cost = get_angle_cost(first_slope, second_slope)[1]
                    second_slope_key = key + f"_{second_slope[0]}_{second_slope[1]}_"
                    second_node_pos_key = nodes_to_nums[second_slope_key + "1"]
                    second_node_neg_key = nodes_to_nums[second_slope_key + "-1"]

                    # Add edges
                    graph.add_edge(second_node_pos_key, first_node_pos_key, other_angle_cost)
                    graph.add_edge(first_node_neg_key, second_node_neg_key, other_angle_cost)

            # Second loop
            for i in range(middle_slopes_length - 1):
                first_slope = slopes["middle"][i]
                first_slope_key = key + f"_{first_slope[0]}_{first_slope[1]}_"
                first_node_pos_key = nodes_to_nums[first_slope_key + "1"]
                first_node_neg_key = nodes_to_nums[first_slope_key + "-1"]
                for j in range(i + 1, middle_slopes_length):
                    second_slope = slopes["middle"][j]
                    costs = get_angle_cost(first_slope, second_slope)
                    angle_cost = costs[0]
                    other_angle_cost = costs[1]
                    second_slope_key = key + f"_{second_slope[0]}_{second_slope[1]}_"
                    second_node_pos_key = nodes_to_nums[second_slope_key + "1"]
                    second_node_neg_key = nodes_to_nums[second_slope_key + "-1"]

                    # Add edges
                    graph.add_edge(first_node_pos_key, second_node_neg_key, angle_cost)
                    graph.add_edge(second_node_neg_key, first_node_pos_key, angle_cost)
                    graph.add_edge(first_node_neg_key, second_node_pos_key, angle_cost)
                    graph.add_edge(second_node_pos_key, first_node_neg_key, angle_cost)

                    graph.add_edge(first_node_pos_key, second_node_pos_key, other_angle_cost)
                    graph.add_edge(second_node_pos_key, first_node_pos_key, other_angle_cost)
                    graph.add_edge(first_node_neg_key, second_node_neg_key, other_angle_cost)
                    graph.add_edge(second_node_neg_key, first_node_neg_key, other_angle_cost)

                for second_slope in end_slopes:
                    costs = get_angle_cost(first_slope, second_slope)
                    angle_cost = costs[0]
                    other_angle_cost = costs[1]            
                    second_slope_key = key + f"_{second_slope[0]}_{second_slope[1]}_"
                    second_node_pos_key = nodes_to_nums[second_slope_key + "1"]
                    second_node_neg_key = nodes_to_nums[second_slope_key + "-1"]

                    # Add edges
                    graph.add_edge(first_node_pos_key, second_node_neg_key, angle_cost)
                    graph.add_edge(second_node_pos_key, first_node_neg_key, angle_cost)

                    graph.add_edge(first_node_neg_key, second_node_neg_key, other_angle_cost)
                    graph.add_edge(second_node_pos_key, first_node_pos_key, other_angle_cost)

            # Third loop
            for i in range(end_slopes_length - 1):
                first_slope = slopes["end"][i]
                first_slope_key = key + f"_{first_slope[0]}_{first_slope[1]}_"
                first_node_pos_key = nodes_to_nums[first_slope_key + "1"]
                first_node_neg_key = nodes_to_nums[first_slope_key + "-1"]
                for j in range(i + 1, end_slopes_length):
                    second_slope = slopes["end"][j]
                    angle_cost = get_angle_cost(first_slope, second_slope)[0]
                    second_slope_key = key + f"_{second_slope[0]}_{second_slope[1]}_"
                    second_node_pos_key = nodes_to_nums[second_slope_key + "1"]
                    second_node_neg_key = nodes_to_nums[second_slope_key + "-1"]            

                    # Add edges
                    graph.add_edge(first_node_pos_key, second_node_neg_key, angle_cost)
                    graph.add_edge(second_node_pos_key, first_node_neg_key, angle_cost)
        
        # Solve graph, print path, and visualize (optional)
        try:
            path_info = find_path(graph, first_start_pos_key, first_target_pos_key)
            num_path = path_info[0]
            optimal_value = path_info[3]
        except:
            print("No path found")
            return None
        path_list = []
        num = num_path.pop(0)
        node = nums_to_nodes[num]
        node_chars = node.split("_")
        node_coord = [int(node_chars[0]), int(node_chars[1])]
        path_list.append(node_coord)
        for num in num_path:
            node = nums_to_nodes[num]
            node_chars = node.split("_")
            node_coord = [int(node_chars[0]), int(node_chars[1])]
            if node_coord != path_list[-1]:
                path_list.append(node_coord)

        print(f"Optimal path: {path_list}")
        print(f"Optimal value: {optimal_value}")
        runtime = time.time() - start_time
        print(f"runtime (s): {runtime}")

        # Visualize
        if visualize == "yes":
            row, column = size, size
            fig, ax = plt.subplots(1, 1)
            x_coords = [point[0] for point in path_list]
            y_coords = [point[1] for point in path_list]
            x = np.linspace(1, row, row)
            y = np.linspace(1, column, column)
            x_grid, y_grid = np.meshgrid(x, y)

            ax.plot(x_grid, y_grid, marker='o', color='k', linestyle='none', markersize=0.2)
            ax.plot(start[0], start[1], marker='o', color='b', linestyle='none', markersize=1)
            ax.plot(target[0], target[1], marker='o', color='r', linestyle='none', markersize=1)
            for edge in edges:
                ax.plot(
                    [edge[0][0], edge[1][0]],
                    [edge[0][1], edge[1][1]],
                    color='k',
                    linewidth=0.2
                )
            ax.plot(x_coords, y_coords, color="g", linewidth=0.2)
            ax.axis("scaled")
            plt.show()

"""
Tổng kết: Quá trình thực hiện cho model 2
Model 2 sử dụng hàm tính năng thêm cost function có sẵn trong thư viện
dijkstar. Đây là hàm phụ thuộc vào một cạnh trước đó và cạnh hiện tại,
chính xác là những gì bài toán này cần. Do đó, đồ thị cần sinh ra ở phần
này nhỏ hơn rất nhiều so với model 1. Các bước tiến hành cũng tương tự:
1. Đọc thông tin của lưới ô vuông
2. (Tuỳ chọn) Xoá đi các đỉnh là đầu mút của tường. Chú ý thiết kế code
sao cho việc xoá bằng list.remove() không trả về lỗi khi không có phần
tử cần xoá.
3. Tạo ra dữ liệu về nodes_info (dictionary) (gồm các cặp đỉnh và hệ số góc tối giản
có thể đi).
4. Từ dữ liệu nodes_info trên, tiếp tục xoá đi các đỉnh "không quan trọng",
là các đỉnh mà chẳng có một hướng nào đi được.
5. Từ active_nodes_info, tạo ra các đỉnh của đồ thị và một quan hệ thứ tự
để phiên dịch từ các đỉnh trừu tượng ra con số và ngược lại. (Quan hệ
này đơn giản hơn ở model 1.)
6. Thêm các cạnh. Chỉ có một kiểu cạnh, là cạnh nối giữa hai đỉnh kề nhau.
7. Giải đồ thị, in ra kết quả về các đỉnh. Ở đây không có sự trùng lặp về đỉnh.
Có thể tích hợp luôn phần vẽ kết quả cuối cùng, hoặc lưu ra một file rồi
vẽ sau.
"""

def solve_with_second_model(size, index, visualize = "no"):
    start_time = time.time()
    # List of all points (to be filtered and used later)
    # Note: j before i
    points_list = [[i, j] for j in range(1, size + 1) for i in range(1, size + 1)]

    nodes_info = {
        f"{i}_{j}": {
            "begin": [], 
            "middle": [], 
            "end": [],
            "reach": [],
        } 
        for j in range(1, size + 1) for i in range(1, size + 1)
    }

    # Two options for nodes_slopes, one through local save, another through
    # in-program generation. Testing seems to prefer the save and read option.
    slopes_path = Path(__file__).parent/"Slopes"/f"Size{size}"
    with open(slopes_path, "r") as f:
        node_slopes = json.load(f)

    # node_slopes = generate_slopes(size=size)

    grid_path = Path(__file__).parent/"Samples"/f"Size{size}"/f"sample{index}.json"
    with open(grid_path, "r") as f:
        grid_info = json.load(f)

    # Get edges, then sort (for faster execution)
    edges = grid_info["edges"]
    edges = sort_edge_list(edge_list=edges)

    # Filter redundant points from points_list
    for edge in edges:
        points_to_remove = get_points_between(edge[0], edge[1])
        for point in points_to_remove:
            # remove from points_list
            try:
                points_list.remove(point)
            except:
                pass

            # remove from nodes_info
            key = f"{point[0]}_{point[1]}"
            try:
                del nodes_info[key]
            except:
                pass
    
    # First validity check
    start = grid_info["start"]
    target = grid_info["target"]
    if start not in points_list: 
        print(f"Invalid: Start is on the wall")
        return None
    if target not in points_list: 
        print(f"Invalid: Target is on the wall")
        return None

    # Nodes_info
    for current_col, current_row in points_list:
        point = [current_col, current_row]
        for slope in node_slopes[f"{current_col}_{current_row}"]:
            max_reach = get_furthest_reach(point, slope, size, size, edges)
            if max_reach > 0:
                nodes_info[f"{current_col}_{current_row}"]["begin"].append(slope)
                nodes_info[f"{current_col}_{current_row}"]["reach"].append(max_reach)
                for num in range(1, max_reach):
                    new_col = current_col + num * slope[0]
                    new_row = current_row + num * slope[1]
                    try:
                        node_slopes[f"{new_col}_{new_row}"].remove(slope)
                    except: 
                        pass
                    nodes_info[f"{new_col}_{new_row}"]["middle"].append(slope)
                new_col = current_col + max_reach * slope[0]
                new_row = current_row + max_reach * slope[1]
                try:
                    node_slopes[f"{new_col}_{new_row}"].remove(slope)
                except:
                    pass
                nodes_info[f"{new_col}_{new_row}"]["end"].append(slope)

    # Second validity check (start and target both in nodes_info)
    if is_empty_dict(nodes_info[f"{start[0]}_{start[1]}"]):
        print(f"Invalid: Start is disconnected")
        return None

    if is_empty_dict(nodes_info[f"{target[0]}_{target[1]}"]):
        print(f"Invalid: Target is disconnected")
        return None
    
    else:
        # Create graph
        graph = Graph()

        # Further filter unreachable nodes from list
        for current_col, current_row in points_list:
            key = f"{current_col}_{current_row}"
            if is_empty_dict(nodes_info[key]):
                del nodes_info[key]
        
        # Crate adjacent edges
        for key in nodes_info.keys():
            point = split_dict_key(key)
            length = len(nodes_info[key]["begin"])
            for i in range(length):
                slope = nodes_info[key]["begin"][i]
                reach = nodes_info[key]["reach"][i]
                unit_length = get_Cartesian_length(slope)
                for j in range(reach):
                    first_point = coords_to_num(point[0] + j * slope[0], point[1] + j * slope[1], size)
                    second_point = coords_to_num(point[0] + (j + 1) * slope[0], point[1] + (j + 1) * slope[1], size)

                    # Add edges
                    graph.add_edge(first_point, second_point, [unit_length, slope])
                    graph.add_edge(second_point, first_point, [unit_length, [-slope[0], -slope[1]]])

        # Solve graph
        start_num = coords_to_num(*start, size)
        target_num = coords_to_num(*target, size)
        try:
            path_info = find_path(graph, start_num, target_num, cost_func=cost_function)
            nums_path = path_info[0]
            optimal_value = path_info[3]
        except:
            print("No path found")
            return None
        points_path = []
        for num in nums_path:
            points_path.append(num_to_coords(num, size))
        print(f"Optimal path: {points_path}")
        print(f"Optimal value {optimal_value}")
        runtime = time.time() - start_time
        print(f"runtime (s): {runtime}")

        # Visualize
        if visualize == "yes":
            row, column = size, size
            fig, ax = plt.subplots(1, 1)
            x_coords = [point[0] for point in points_path]
            y_coords = [point[1] for point in points_path]
            x = np.linspace(1, row, row)
            y = np.linspace(1, column, column)
            x_grid, y_grid = np.meshgrid(x, y)

            ax.plot(x_grid, y_grid, marker='o', color='k', linestyle='none', markersize=0.2)
            ax.plot(start[0], start[1], marker='o', color='b', linestyle='none', markersize=1)
            ax.plot(target[0], target[1], marker='o', color='r', linestyle='none', markersize=1)
            for edge in edges:
                ax.plot(
                    [edge[0][0], edge[1][0]],
                    [edge[0][1], edge[1][1]],
                    color='k',
                    linewidth=0.2
                )
            ax.plot(x_coords, y_coords, color="g", linewidth=0.2)
            ax.axis('scaled')
            plt.show()