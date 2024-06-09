from pathlib import Path
import json
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import math
import time

class Constant:
    """
    Turning time: After normalizing grid size to 1 and velcity to 1.
    Additional assumption: Turning time for acute angles 
    is twice as long as that for obtuse angles.
    """
    # Based on information in the video "The Fastest Maze-Solving Competition On Earth"
    # of channel Veritasium, from 10:40 to 11:20.
    acute = 7.5
    obtuse_or_right = 3.75

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

# Get points on a segment, for first-stage filtering of redundant vertices,
# among other purposes
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

# Get maximum number of vertices for a feasible path of given maze
def get_maximum_number_of_vertices(size, index):
    """
    Get maximum number of vertices for a feasible path of given maze.
    """
    grid_path = Path(__file__).parent/"Samples"/f"Size{size}"/f"sample{index}.json"
    with open(grid_path, "r") as f:
        edges = json.load(f)["edges"]
    redundant_points = set()
    # Filter redundant points from points_list
    for edge in edges:
        points_to_remove = get_points_between(edge[0], edge[1])
        for point in points_to_remove:
            # add to redundant_points
            redundant_points.add(tuple(point))
    redundant_points = list(redundant_points)
    return size**2 - len(redundant_points)

# Sort edge_list, for faster execution
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

# Generate unreachable_nodes
def get_unreachable_nodes(size, index):
    """
    In-program generation of unreachable_nodes
    """
    # List of all points (to be filtered and used later)
    # Note: j before i
    points_list = [[i, j] for j in range(1, size + 1) for i in range(1, size + 1)]

    unreachable_nodes = {
        f"{i}_{j}": []
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
    redundant_points = set()
    # Filter redundant points from points_list
    for edge in edges:
        points_to_remove = get_points_between(edge[0], edge[1])
        for point in points_to_remove:
            # remove from points_list
            try:
                points_list.remove(point)
            except:
                pass

            # remove from unreachable_nodes
            key = f"{point[0]}_{point[1]}"
            try:
                del unreachable_nodes[key]
            except:
                pass
            # add to redundant_points
            redundant_points.add(tuple(point))
    redundant_points = list(redundant_points)

    # Generate unreachable_nodes
    for current_col, current_row in points_list:
        point = [current_col, current_row]
        for slope in node_slopes[f"{current_col}_{current_row}"]:
            # The following is to get the furthest possible reach in a
            # given direction
            col_slope, row_slope = slope[0], slope[1]
            # Get maximum column reach
            if col_slope == 0: max_col_reach = size
            elif col_slope > 0: max_col_reach = int((size - current_col) / col_slope) + 1
            else: max_col_reach = int((current_col - 1) / (-col_slope)) + 1
            # Get maximum row reach
            if row_slope == 0: max_row_reach = size
            else: max_row_reach = int((size - current_row) / row_slope) + 1

            min_reach = 0
            max_reach = min(max_row_reach, max_col_reach)
            max_valid_reach = max_reach - 1
            while max_reach - min_reach > 1:
                temp = min_reach
                new_reach = int((max_reach + min_reach) / 2)
                min_reach = new_reach
                test_point = [current_col + new_reach * col_slope, current_row + new_reach * row_slope]
                for edge in edges:
                    if check_intersection(point, test_point, edge):
                        max_reach = new_reach
                        min_reach = temp
                        break

            first_group = []
            second_group = []

            for i in range(min_reach + 1):
                point_to_add = [current_col + i * col_slope, current_row + i * row_slope]
                if tuple(point_to_add) not in redundant_points:
                    first_group.append(point_to_add)

            for i in range(min_reach + 1, max_valid_reach + 1):
                point_to_add = [current_col + i * col_slope, current_row + i * row_slope]
                if tuple(point_to_add) not in redundant_points:
                    second_group.append(point_to_add)

            for checkpoint in second_group:
                key = f"{checkpoint[0]}_{checkpoint[1]}"
                unreachable_nodes[key] += first_group
                
            first_point = first_group.pop(0)
            key = f"{first_point[0]}_{first_point[1]}"
            unreachable_nodes[key] += second_group
            for checkpoint in first_group:
                key = f"{checkpoint[0]}_{checkpoint[1]}"
                unreachable_nodes[key] += second_group
                try:
                    node_slopes[key].remove(slope)
                except:
                    pass
    return unreachable_nodes

# Hàm giải với số bước cố định chọn trước
def solve_maze_with_given_step(size, index, step, status = "optimal"):
    """
    Solve given maze for a solution with fixed number of steps, optimal
    if 'status' = "optimal", else an arbitrary feasible solution.
    """
    # Thiết lập mô hình
    model = gp.Model()

    # Look for any feasible solution if 'status' = "feasible"
    if status == "feasible": model.params.SolutionLimit = 1

    # Generation of nodes through pre-written data, only available for
    # existing samples of size 40*40 and below. In the long run, an
    # in-program generation is prefered.
    nodes_path = Path(__file__).parent/"Unreachable_nodes"/f"Size{size}"/f"sample{index}.json"
    with open(nodes_path, "r") as f:
        unreachable_nodes = json.load(f)

    # In-program generation
    # unreachable_nodes = get_unreachable_nodes(size, index)

    grid_path = Path(__file__).parent/"Samples"/f"Size{size}"/f"sample{index}.json"
    with open(grid_path, "r") as f:
        maze = json.load(f)
    start = maze["start"]
    target = maze["target"]

    # Get edges
    edges = maze["edges"]

    redundant_points = set()
    points_list = [[i, j] for i in range(1, size + 1) for j in range(1, size + 1)]
    # Filter redundant points from points_list
    for edge in edges:
        points_to_remove = get_points_between(edge[0], edge[1])
        for point in points_to_remove:
            # remove from points_list
            try:
                points_list.remove(point)
            except:
                pass

            # add to redundant_points
            redundant_points.add(tuple(point))
    redundant_points = list(redundant_points)
    
    # Check validity of number of steps
    max_step = size**2 - len(redundant_points)
    if step > max_step: 
        raise ValueError(f"Số bước vượt quá số bước tối đa: {max_step}")

    n = size
    M = 2*n**2 + 10
    N = step

    x = np.empty((n+1, n+1, N+1), dtype=object)
    for i in range(1,n+1):
        for j in range(1,n+1):
            for k in range(1, N+1):
                x[i,j,k] = model.addVar(vtype=GRB.BINARY)

    # Tất cả các đỉnh trong redundant_points đều bằng 0
    for k in range(1, N + 1):
        for i, j in redundant_points:
            model.addConstr(x[i, j, k] == 0)

    # Tính tọa độ a, b của đỉnh thứ k
    a = np.empty((N+1), dtype=object)
    b = np.empty((N+1), dtype=object)
    for k in range(1, N+1):
        a[k] = model.addVar(lb=0, ub=n, vtype= GRB.CONTINUOUS)
        model.addConstr(gp.quicksum([x[i, j, k]*i for i in range(1,n+1) for j in range(1,n+1)]) == a[k])
        b[k] = model.addVar(lb=0, ub=n, vtype= GRB.CONTINUOUS)
        model.addConstr(gp.quicksum([x[i, j, k]*j for i in range(1,n+1) for j in range(1,n+1)]) == b[k])

    # Biến kiểm tra xem đỉnh thứ k có được chọn hay không
    count = np.empty(N + 1, dtype=object)
    for k in range(1, N + 1):
        count[k] = model.addVar(vtype=GRB.BINARY)
        model.addConstr(count[k] == gp.quicksum([x[i, j, k] for i, j in points_list]))

    # Tính quãng đường
    # Không nhân được ba lần
    step_length_squared = np.empty(N, dtype=object)
    true_step_length_squared = np.empty(N, dtype=object)
    sqrt_var = np.empty(N, dtype=object)
    for k in range(1, N):
        step_length_squared[k] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.INTEGER)
        model.addConstr(step_length_squared[k] == (a[k+1] - a[k])**2 + (b[k+1] - b[k])**2)
        true_step_length_squared[k] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.INTEGER)
        model.addConstr(true_step_length_squared[k] == step_length_squared[k] * count[k + 1])
        sqrt_var[k] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="sqrt_var")
        model.addConstr(sqrt_var[k] * sqrt_var[k] == true_step_length_squared[k], "sqrt_constr")
    sqrt = gp.quicksum([sqrt_var[k] for k in range(1, N)])

    # Danh sách các véc-tơ giữa hai đỉnh liên tiếp
    vectors_list = []
    for k in range(1, N):
        vectors_list.append([a[k + 1] - a[k], b[k + 1] - b[k]])

    # Dùng luôn từ vectors_list, tránh tính lại nhiều lần trong hàm dot_product.
    # Hơn nữa vectors_list còn được dùng trong điều kiện không thẳng hàng ở sau.
    total_angle_cost = 0
    for i in range(2, N):
        vector1, vector2 = vectors_list[i - 2], vectors_list[i - 1]
        dot_product = vector1[0]*vector2[0] + vector1[1]*vector2[1]
        b1 = model.addVar(vtype = GRB.BINARY, name = "b1") # = 1 nếu dot_product > 0, = 0 nếu ngược lại
        b2 = model.addVar(vtype = GRB.BINARY, name = "b2") # = 1 nếu dot_product < 0, = 0 nếu ngược lại
        b3 = model.addVar(vtype = GRB.BINARY, name = "b3") # = 1 nếu dot_product = 0

        # Chỉ chính xác hoá b1, b2. Điều kiện sau sẽ làm chính xác b3.
        model.addConstr(b1+b2+b3 == 1)
        model.addConstr(dot_product <= M*b1)
        model.addConstr(dot_product >= (M + 1)*(b1 - 1) + 1)
        model.addConstr(dot_product >= -M*b2)
        model.addConstr(dot_product <= (M + 1)*(1 - b2) - 1)

        # Tính thời gian quay, cộng vào tổng
        angle_cost = (Constant.obtuse_or_right * (b1 + b3) + Constant.acute * b2) * count[i + 1]
        total_angle_cost += angle_cost

    # Hàm mục tiêu
    obj = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
    model.addConstr(obj == sqrt + total_angle_cost)
    model.setObjective(obj, GRB.MINIMIZE)
    # model.setObjective(0, sense=GRB.MINIMIZE)

    # Các ràng buộc
    # Mỗi bước thứ k, chọn đúng một đỉnh
    for k in range(1, N + 1):
        model.addConstr(gp.quicksum([x[i, j, k] for i, j in points_list]) == 1)

    # Mỗi đỉnh (i, j) xuất hiện nhiều nhất 1 lần
    for i, j in points_list:
        model.addConstr(gp.quicksum([x[i, j, k] for k in range(1, N + 1)]) <= 1)

    """
    Vấn đề tiếp theo là tính góc trong trường hợp có ba đỉnh thẳng hàng, ở
    đây góc sẽ bằng 0, hoặc 180 độ.
    Hướng 1: Thêm điều kiện để loại bỏ ba đỉnh thẳng hàng.
    Hướng 2: Xử lý góc riêng biệt (chú ý trường hợp 180 độ)???
    """

    # Hướng 1: Không cho phép ba đỉnh liên tiếp thẳng hàng
    # link: https://or.stackexchange.com/questions/7726/no-not-equals-constraint-in-gurobi
    """
    If you want x1 ≠ x2, you can linearize |x1 - x2| ≥ 𝜀, where 𝜀
    is your tolerance.
    You can do this by introducing a boolean variable y and by imposing:
    x1 - x2 ≤ -𝜀 + My and x1 - x2 ≥ 𝜀 - (1 - y)M,
    where M is the smallest large constant you can think of.
    For integer formulations, 𝜀 = 1 and M = (max of x1) - (min of x2) + 1
    """
    for i in range(N - 2):
        check_collinear = vectors_list[i][0] * vectors_list[i + 1][1] - vectors_list[i][1] * vectors_list[i + 1][0]
        # Điều kiện: 'check_collinear' khác 0. Tham khảo link 
        # (lưu ý: check_collinear là một số nguyên)
        b = model.addVar(vtype=GRB.BINARY)
        model.addConstr(check_collinear <= -1 + M*b)
        model.addConstr(check_collinear >= 1 - M*(1 - b))

    # Toạ độ đỉnh đầu, cuối
    model.addConstr(x[start[0],start[1],1] == 1)
    model.addConstr(x[target[0],target[1],N] == 1)

    # Điều kiện không chạm tường, được xử lý trước bằng một bước riêng biệt
    for i, j in points_list:
        unreachable_points = unreachable_nodes[f"{i}_{j}"]
        for k in range(1, N):
            # Mô tả điều kiện: Nếu x[i, j, k] = 1 thì tổng các x[i', j', k + 1] bằng 0,
            # với (i', j') là một đỉnh trong 'unreachable_points' của (i, j)
            model.addConstr(x[i, j, k] + gp.quicksum([x[point[0], point[1], k + 1] for point in unreachable_points]) <= 1)

    # Các điều kiện sau là các điều kiện luôn đúng, nhưng chúng có ảnh hưởng
    # đến quá trình giải (theo cách tốt hoặc xấu). Vì thế, hãy thử thêm hoặc
    # bỏ các điều kiện này trong quá trình giải.

    # Có ít nhất 2 bước
    # model.addConstr(count[2] == 1)

    # Chặn trên đã biết cho hàm mục tiêu (nếu có)
    # model.addConstr(obj <= )

    # Tìm nghiệm tối ưu
    model.optimize()

    # In ra thông tin về quãng đường, dùng cho bước hai
    if model.status == GRB.OPTIMAL or model.status == GRB.SOLUTION_LIMIT:
        return True, model.ObjVal
    else: return False, 0

# Hàm giải với số bước không quá một chặn trên cho trước
def solve_for_solution_with_bounded_steps(size, index, step_bound, status = "optimal"):
    """
    Solve given maze for solution, optimal if 'status' = "optimal", else
    an arbitrary feasible solution.
    """
    # Thiết lập mô hình
    model = gp.Model()

    # Look for any feasible solution if 'status' = "feasible"
    if status == "feasible": model.params.SolutionLimit = 1

    # Cần để làm phép nhân ba biến
    model.params.NonConvex = 2

    # Lấy thông tin
    nodes_path = Path(__file__).parent/"Unreachable_nodes"/f"Size{size}"/f"sample{index}.json"
    grid_path = Path(__file__).parent/"Samples"/f"Size{size}"/f"sample{index}.json"
    with open(nodes_path, "r") as f:
        unreachable_nodes = json.load(f)

    # In-program generation
    # unreachable_nodes = get_unreachable_nodes(size, index)

    with open(grid_path, "r") as f:
        maze = json.load(f)
    start = maze["start"]
    target = maze["target"]

    # Lấy tập cạnh là tường
    edges = maze["edges"]

    redundant_points = set()
    points_list = [[i, j] for i in range(1, size + 1) for j in range(1, size + 1)]
    # Loại những đỉnh thừa
    for edge in edges:
        points_to_remove = get_points_between(edge[0], edge[1])
        for point in points_to_remove:
            try:
                points_list.remove(point)
            except:
                pass

            # Thêm đỉnh thừa vào một danh sách
            redundant_points.add(tuple(point))
    redundant_points = list(redundant_points)

    # Check validity of number of steps
    max_step = size**2 - len(redundant_points)
    if step_bound > max_step: 
        raise ValueError(f"Số bước vượt quá số bước tối đa: {max_step}")

    n = size # cho gọn trong mô hình

    # Hằng số chặn trên, dùng cho các ràng buộc "nếu, thì"
    M = 2*n**2 + 10

    # Chặn trên cho số bước tối đa
    N = step_bound

    # Lập mô hình
    model = gp.Model()
    x = np.empty((n+1, n+1, N+1), dtype=object)
    for i in range(1,n+1):
        for j in range(1,n+1):
            for k in range(1, N+1):
                x[i,j,k] = model.addVar(vtype=GRB.BINARY)

    # Tất cả các đỉnh trong redundant_points đều bằng 0
    for k in range(1, N + 1):
        for i, j in redundant_points:
            model.addConstr(x[i, j, k] == 0)

    # Tính tọa độ a, b của đỉnh thứ k
    a = np.empty((N+1), dtype=object)
    b = np.empty((N+1), dtype=object)
    for k in range(1, N+1):
        a[k] = model.addVar(lb=0, ub=n, vtype= GRB.CONTINUOUS)
        model.addConstr(gp.quicksum([x[i, j, k]*i for i in range(1,n+1) for j in range(1,n+1)]) == a[k])
        b[k] = model.addVar(lb=0, ub=n, vtype= GRB.CONTINUOUS)
        model.addConstr(gp.quicksum([x[i, j, k]*j for i in range(1,n+1) for j in range(1,n+1)]) == b[k])

    # Biến kiểm tra xem đỉnh thứ k có được chọn hay không
    count = np.empty(N + 1, dtype=object)
    for k in range(1, N + 1):
        count[k] = model.addVar(vtype=GRB.BINARY)
        model.addConstr(count[k] == gp.quicksum([x[i, j, k] for i, j in points_list]))

    # Tính quãng đường
    # Không nhân được ba lần
    step_length_squared = np.empty(N, dtype=object)
    true_step_length_squared = np.empty(N, dtype=object)
    sqrt_var = np.empty(N, dtype=object)
    for k in range(1, N):
        step_length_squared[k] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.INTEGER)
        model.addConstr(step_length_squared[k] == (a[k+1] - a[k])**2 + (b[k+1] - b[k])**2)
        true_step_length_squared[k] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.INTEGER)
        model.addConstr(true_step_length_squared[k] == step_length_squared[k] * count[k + 1])
        sqrt_var[k] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="sqrt_var")
        model.addConstr(sqrt_var[k] * sqrt_var[k] == true_step_length_squared[k], "sqrt_constr")
    sqrt = gp.quicksum([sqrt_var[k] for k in range(1, N)])

    # Danh sách các véc-tơ giữa hai đỉnh liên tiếp
    vectors_list = []
    for k in range(1, N):
        vectors_list.append([a[k + 1] - a[k], b[k + 1] - b[k]])

    # Tính tích vô hướng, qua đó xác định giá trị của góc quay.
    # Dùng luôn từ vectors_list, tránh tính lại nhiều lần trong hàm dot_product.
    # Hơn nữa vectors_list còn được dùng trong điều kiện không thẳng hàng ở sau.
    total_angle_cost = 0
    for i in range(2, N):
        vector1, vector2 = vectors_list[i - 2], vectors_list[i - 1]
        dot_product = vector1[0]*vector2[0] + vector1[1]*vector2[1]
        # angle[i] = model.addVar(lb=Constant.obtuse_or_right, ub=Constant.acute, vtype=GRB.CONTINUOUS)
        b1 = model.addVar(vtype = GRB.BINARY, name = "b1") # = 1 nếu dot_product > 0, = 0 nếu ngược lại
        b2 = model.addVar(vtype = GRB.BINARY, name = "b2") # = 1 nếu dot_product < 0, = 0 nếu ngược lại
        b3 = model.addVar(vtype = GRB.BINARY, name = "b3") # = 1 nếu dot_product = 0

        # Chỉ chính xác hoá b1, b2. Điều kiện sau sẽ làm chính xác b3.
        model.addConstr(b1+b2+b3 == 1)
        model.addConstr(dot_product <= M*b1)
        model.addConstr(dot_product >= (M + 1)*(b1 - 1) + 1)
        model.addConstr(dot_product >= -M*b2)
        model.addConstr(dot_product <= (M + 1)*(1 - b2) - 1)

        # Tính thời gian quay, cộng vào tổng
        angle_cost = (Constant.obtuse_or_right * (b1 + b3) + Constant.acute * b2) * count[i + 1]
        total_angle_cost += angle_cost

    # Hàm mục tiêu
    obj = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
    model.addConstr(obj == sqrt + total_angle_cost)
    model.setObjective(obj, GRB.MINIMIZE)
    # model.setObjective(0, sense=GRB.MINIMIZE)

    # Các ràng buộc
    # Mỗi bước thứ k, chọn nhiều nhất một đỉnh
    for k in range(1, N + 1):
        model.addConstr(count[k] <= 1)

    # Mỗi đỉnh (i, j) xuất hiện nhiều nhất một lần
    for i, j in points_list:
        model.addConstr(gp.quicksum([x[i, j, k] for k in range(1, N + 1)]) <= 1)

    # Nếu đỉnh thứ k được chọn thì tất cả các đỉnh trước đó cũng được chọn
    for k in range(1, N):
        model.addConstr(count[k] >= count[k + 1])

    for i in range(N - 2):
        check_collinear = vectors_list[i][0] * vectors_list[i + 1][1] - vectors_list[i][1] * vectors_list[i + 1][0]
        # Điều kiện: 'check_collinear' khác 0. Tham khảo link https://or.stackexchange.com/questions/7726/no-not-equals-constraint-in-gurobi
        # (lưu ý: check_collinear là một số nguyên)
        b = model.addVar(vtype=GRB.BINARY)
        # Lưu ý: Quan tâm đến đỉnh thứ i + 3. Nếu nó không được chọn, điều kiện này luôn được thoả mãn với b = 0
        model.addConstr(check_collinear - M*(1 - count[i + 3]) <= -1 + M*b)
        model.addConstr(check_collinear >= 1 - M*(1 - b))

    # Toạ độ đỉnh đầu, cuối
    # Lưu ý về đỉnh cuối!
    model.addConstr(x[start[0],start[1],1] == 1)  
    model.addConstr(gp.quicksum([x[target[0], target[1], k] for k in range(1, N + 1)]) == 1)

    # Vị trí của đỉnh cuối phải thực sự là vị trí cuối
    model.addConstr(gp.quicksum([x[target[0], target[1], k]*k for k in range(1, N + 1)]) == gp.quicksum([count[k] for k in range(1, N + 1)]))

    # Điều kiện không chạm tường, được xử lý trước bằng một bước riêng biệt
    for i, j in points_list:
        unreachable_points = unreachable_nodes[f"{i}_{j}"]
        for k in range(1, N):
            # Mô tả điều kiện: Nếu x[i, j, k] = 1 thì tổng các x[i', j', k + 1] bằng 0,
            # với (i', j') là một đỉnh trong 'unreachable_points' của (i, j)
            model.addConstr(x[i, j, k] + gp.quicksum([x[point[0], point[1], k + 1] for point in unreachable_points]) <= 1)

    # Các điều kiện sau là các điều kiện luôn đúng, nhưng chúng có ảnh hưởng
    # đến quá trình giải (theo cách tốt hoặc xấu). Vì thế, hãy thử thêm hoặc
    # bỏ các điều kiện này trong quá trình giải.

    # Có ít nhất 2 bước
    # model.addConstr(count[2] == 1)

    # Chặn trên đã biết cho hàm mục tiêu (nếu có)
    # model.addConstr(obj <= )

    model.optimize()

    # In ra thông tin về quãng đường, dùng cho bước hai
    if model.status == GRB.OPTIMAL or model.status == GRB.SOLUTION_LIMIT:
        return True, model.ObjVal
    else: return False, 0   

# Hàm giải "tổng quát"
def solve_maze(size, index, bound_for_feasibility = None, status_for_feasibility = "optimal", method: int = 2):
    """
    Solve given maze.
    NOTE:
    - If bound_for_feasibility = None, model will be constructed with
    default upper bound for maximum number of vertices, which may result
    in a model with O(size**4) constraints, and this is too much for
    size > 10.
    - Else, the process will start searching for a feasible solution with
    an upper bound of vertices given by 'bound_for_feasibility', using the 
    fixed-step version ('method' = 1) or bounded step version ('method' = 2).
    If no feasible solution is found, the process is stopped, otherwise
    there will be a much better bound for maximum number of vertices
    obtained from the feasible objective value found above, thereby 
    decreasing the size of the final model, raising solvability.
    """
    start_time = time.time()
    # Thiết lập mô hình
    model = gp.Model()

    # Cần để làm phép nhân ba biến
    model.params.NonConvex = 2

    # Lấy thông tin
    nodes_path = Path(__file__).parent/"Unreachable_nodes"/f"Size{size}"/f"sample{index}.json"
    grid_path = Path(__file__).parent/"Samples"/f"Size{size}"/f"sample{index}.json"
    with open(nodes_path, "r") as f:
        unreachable_nodes = json.load(f)

    # In-program generation
    # unreachable_nodes = get_unreachable_nodes(size, index)

    with open(grid_path, "r") as f:
        maze = json.load(f)
    start = maze["start"]
    target = maze["target"]

    # Lấy tập cạnh là tường
    edges = maze["edges"]

    redundant_points = set()
    points_list = [[i, j] for i in range(1, size + 1) for j in range(1, size + 1)]
    # Loại những đỉnh thừa
    for edge in edges:
        points_to_remove = get_points_between(edge[0], edge[1])
        for point in points_to_remove:
            try:
                points_list.remove(point)
            except:
                pass

            # Thêm đỉnh thừa vào một danh sách
            redundant_points.add(tuple(point))
    redundant_points = list(redundant_points)

    # Đầu tiên, tìm một nghiệm tối ưu cho số bước cụ thể, là một nghiệm chấp nhận được
    # cho bài toán với số bước chưa biết. Từ giá trị hàm mục tiêu ở đây, ta
    # thu được một chặn trên cho số bước tối đa.
    max_step = size**2 - len(redundant_points)
    if bound_for_feasibility > max_step:
        raise ValueError(f"Số bước vượt quá số bước tối đa: {max_step}")

    elif bound_for_feasibility == None: N = max_step

    # Chặn trên cho số bước tối đa, từ việc độ dài mỗi đoạn ít nhất là 1 đơn
    # vị, và vì không có ba đỉnh liên tiếp thẳng hàng, mỗi bước đi đều bao
    # gồm một lần quay. Do đó, mỗi bước đi sẽ tốn ít nhất 1 + 3.75 = 4.75 đvtg.đvtg
    # Có N - 1 bước đi như vậy.
    else:
        print(f"*\n*\n*\n*\n*")
        print(f"Begin step one: Finding a {status_for_feasibility} solution within range of {bound_for_feasibility} vertices using method {method}")
        if method == 1:
            count = 2
            for step in range(2, bound_for_feasibility + 1):
                info = solve_maze_with_given_step(size, index, step, status_for_feasibility)
                if info[0]:
                    N = int(info[1]/4.75) + 2
                    print(f"*\n*\n*\n*\n*")
                    print(f"Obtained a bound for maximum number of vertices: {N}")
                    print(f"Begin step two: Solving given maze with maximum {N} vertices")
                    break
                count += 1
            if count == bound_for_feasibility + 1:
                print(f"*\n*\n*\n*\n*")
                print(f"No feasible solution found with given range of {bound_for_feasibility} vertices")
                return None
        elif method == 2:
            info = solve_for_solution_with_bounded_steps(size, index, bound_for_feasibility, status_for_feasibility)
            if info[0]:
                N = int(info[1]/4.75) + 2
                print(f"*\n*\n*\n*\n*")
                print(f"Obtained a bound for maximum number of vertices: {N}")
                print(f"Begin step two: Solving given maze with maximum {N} vertices")
            else:
                print(f"*\n*\n*\n*\n*")
                print(f"No feasible solution found with given range of {bound_for_feasibility} vertices")
                return None

    n = size # cho gọn trong mô hình

    # Hằng số chặn trên, dùng cho các ràng buộc "nếu, thì"
    M = 2*n**2 + 10

    # Lập mô hình
    model = gp.Model()
    x = np.empty((n+1, n+1, N+1), dtype=object)
    for i in range(1,n+1):
        for j in range(1,n+1):
            for k in range(1, N+1):
                x[i,j,k] = model.addVar(vtype=GRB.BINARY)

    # Tất cả các đỉnh trong redundant_points đều bằng 0
    for k in range(1, N + 1):
        for i, j in redundant_points:
            model.addConstr(x[i, j, k] == 0)

    # Tính tọa độ a, b của đỉnh thứ k
    a = np.empty((N+1), dtype=object)
    b = np.empty((N+1), dtype=object)
    for k in range(1, N+1):
        a[k] = model.addVar(lb=0, ub=n, vtype= GRB.CONTINUOUS)
        model.addConstr(gp.quicksum([x[i, j, k]*i for i in range(1,n+1) for j in range(1,n+1)]) == a[k])
        b[k] = model.addVar(lb=0, ub=n, vtype= GRB.CONTINUOUS)
        model.addConstr(gp.quicksum([x[i, j, k]*j for i in range(1,n+1) for j in range(1,n+1)]) == b[k])

    # Biến kiểm tra xem đỉnh thứ k có được chọn hay không
    count = np.empty(N + 1, dtype=object)
    for k in range(1, N + 1):
        count[k] = model.addVar(vtype=GRB.BINARY)
        model.addConstr(count[k] == gp.quicksum([x[i, j, k] for i, j in points_list]))

    # Tính quãng đường
    # Không nhân được ba lần
    step_length_squared = np.empty(N, dtype=object)
    true_step_length_squared = np.empty(N, dtype=object)
    sqrt_var = np.empty(N, dtype=object)
    for k in range(1, N):
        step_length_squared[k] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.INTEGER)
        model.addConstr(step_length_squared[k] == (a[k+1] - a[k])**2 + (b[k+1] - b[k])**2)
        true_step_length_squared[k] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.INTEGER)
        model.addConstr(true_step_length_squared[k] == step_length_squared[k] * count[k + 1])
        sqrt_var[k] = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="sqrt_var")
        model.addConstr(sqrt_var[k] * sqrt_var[k] == true_step_length_squared[k], "sqrt_constr")
    sqrt = gp.quicksum([sqrt_var[k] for k in range(1, N)])

    # Danh sách các véc-tơ giữa hai đỉnh liên tiếp
    vectors_list = []
    for k in range(1, N):
        vectors_list.append([a[k + 1] - a[k], b[k + 1] - b[k]])

    # Tính tích vô hướng, qua đó xác định giá trị của góc quay.
    # Dùng luôn từ vectors_list, tránh tính lại nhiều lần trong hàm dot_product.
    # Hơn nữa vectors_list còn được dùng trong điều kiện không thẳng hàng ở sau.
    total_angle_cost = 0
    for i in range(2, N):
        vector1, vector2 = vectors_list[i - 2], vectors_list[i - 1]
        dot_product = vector1[0]*vector2[0] + vector1[1]*vector2[1]
        # angle[i] = model.addVar(lb=Constant.obtuse_or_right, ub=Constant.acute, vtype=GRB.CONTINUOUS)
        b1 = model.addVar(vtype = GRB.BINARY, name = "b1") # = 1 nếu dot_product > 0, = 0 nếu ngược lại
        b2 = model.addVar(vtype = GRB.BINARY, name = "b2") # = 1 nếu dot_product < 0, = 0 nếu ngược lại
        b3 = model.addVar(vtype = GRB.BINARY, name = "b3") # = 1 nếu dot_product = 0

        # Chỉ chính xác hoá b1, b2. Điều kiện sau sẽ làm chính xác b3.
        model.addConstr(b1+b2+b3 == 1)
        model.addConstr(dot_product <= M*b1)
        model.addConstr(dot_product >= (M + 1)*(b1 - 1) + 1)
        model.addConstr(dot_product >= -M*b2)
        model.addConstr(dot_product <= (M + 1)*(1 - b2) - 1)

        # Tính thời gian quay, cộng vào tổng
        angle_cost = (Constant.obtuse_or_right * (b1 + b3) + Constant.acute * b2) * count[i + 1]
        total_angle_cost += angle_cost

    # Hàm mục tiêu
    obj = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
    model.addConstr(obj == sqrt + total_angle_cost)
    model.setObjective(obj, GRB.MINIMIZE)
    # model.setObjective(0, sense=GRB.MINIMIZE)

    # Các ràng buộc
    # Mỗi bước thứ k, chọn nhiều nhất một đỉnh
    for k in range(1, N + 1):
        model.addConstr(count[k] <= 1)

    # Mỗi đỉnh (i, j) xuất hiện nhiều nhất một lần
    for i, j in points_list:
        model.addConstr(gp.quicksum([x[i, j, k] for k in range(1, N + 1)]) <= 1)

    # Nếu đỉnh thứ k được chọn thì tất cả các đỉnh trước đó cũng được chọn
    for k in range(1, N):
        model.addConstr(count[k] >= count[k + 1])

    for i in range(N - 2):
        check_collinear = vectors_list[i][0] * vectors_list[i + 1][1] - vectors_list[i][1] * vectors_list[i + 1][0]
        # Điều kiện: 'check_collinear' khác 0. Tham khảo link https://or.stackexchange.com/questions/7726/no-not-equals-constraint-in-gurobi
        # (lưu ý: check_collinear là một số nguyên)
        b = model.addVar(vtype=GRB.BINARY)
        # Lưu ý: Quan tâm đến đỉnh thứ i + 3. Nếu nó không được chọn, điều kiện này luôn được thoả mãn với b = 0
        model.addConstr(check_collinear - M*(1 - count[i + 3]) <= -1 + M*b)
        model.addConstr(check_collinear >= 1 - M*(1 - b))

    # Toạ độ đỉnh đầu, cuối
    # Lưu ý về đỉnh cuối!
    model.addConstr(x[start[0],start[1],1] == 1)  
    model.addConstr(gp.quicksum([x[target[0], target[1], k] for k in range(1, N + 1)]) == 1)

    # Vị trí của đỉnh cuối phải thực sự là vị trí cuối
    model.addConstr(gp.quicksum([x[target[0], target[1], k]*k for k in range(1, N + 1)]) == gp.quicksum([count[k] for k in range(1, N + 1)]))

    # Điều kiện không chạm tường, được xử lý trước bằng một bước riêng biệt
    for i, j in points_list:
        unreachable_points = unreachable_nodes[f"{i}_{j}"]
        for k in range(1, N):
            # Mô tả điều kiện: Nếu x[i, j, k] = 1 thì tổng các x[i', j', k + 1] bằng 0,
            # với (i', j') là một đỉnh trong 'unreachable_points' của (i, j)
            model.addConstr(x[i, j, k] + gp.quicksum([x[point[0], point[1], k + 1] for point in unreachable_points]) <= 1)

    # Các điều kiện sau là các điều kiện luôn đúng, nhưng chúng có ảnh hưởng
    # đến quá trình giải (theo cách tốt hoặc xấu). Vì thế, hãy thử thêm hoặc
    # bỏ các điều kiện này trong quá trình giải.

    # Có ít nhất 2 bước
    # model.addConstr(count[2] == 1)

    # Chặn trên đã biết cho hàm mục tiêu (nếu có)
    # model.addConstr(obj <= solve_maze_with_given_step(size, index, step)[1]+0.01)

    model.optimize()
    # In ra thông tin về quãng đường
    if model.status == GRB.OPTIMAL:
        print(f"Optimal objective value: {model.objVal}")
        for k in range(1, N + 1):
            if count[k].x == 0: break
            for i, j in points_list:
                if x[i, j, k].x == 1:
                    print(f"{k}_th vertex: ({i}, {j})")
    runtime = time.time() - start_time
    print(f"Runtime (s): {runtime}")
