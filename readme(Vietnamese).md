Xin chào! Đây là bài tập cuối kỳ môn Tối ưu hoá.

## Mô tả bài toán
Bài toán nhóm đặt ra được lấy cảm hứng từ cuộc thi Micromouse (xem trang [Wiki](https://en.wikipedia.org/wiki/Micromouse) này để biết thêm thông tin). Nó được phát biểu như sau:
Cho một lưới đỉnh cỡ m*n với các vật cản ở dạng các cạnh nối hai điểm bất kỳ trên lưới, hãy tìm đường đi nhanh nhất trên lưới nối một cặp điểm (xuất phát-mục tiêu) cho trước. Có hai tiêu chí để đánh giá đường đi như sau
- Độ dài đường đi (tính theo khoảng cách Euclide)
- Số khúc rẽ (chia thành góc rẽ là góc nhọn và ngược lại)

(Vận tốc di chuyển và chi phí rẽ đều được giả định là không đổi). Để biết thêm thông tin về cách xác định các hằng số, các nhận xét và bàn luận về ý tưởng giải bài toán, hãy xem các slide trình bày của nhóm [tại đây]().

Dưới đây là ví dụ về một lưới như vậy, với điểm bắt đầu được đánh dấu màu xanh, mục tiêu màu đỏ và chướng ngại vật màu đen

![](/Sample_images/Size4/sample3.png "Lưới ví dụ")

và đây là một đường đi đề xuất

![](/Proposed_solutions/Visualize/Size4/sample3.png "Đường đi đề xuất")

### Tổng quan về các phương pháp sử dụng
Hai phương pháp đã được nhóm sử dụng để làm bài toán:
- Đầu tiên là thông qua tối ưu phi tuyến, giải bằng gurobi. Nhóm đã xử lý được một phần các mẫu được tạo sẵn, cụ thể là những mẫu có nghiệm chấp nhận được.
- Thứ hai là thông qua đồ thị, giải bằng thuật toán Dijkstra. Nhóm thực hiện trên python với sự trợ giúp từ thư viện [dijkstar](https://pypi.org/project/Dijkstar/).

### Hướng dẫn sử dụng
Các tài nguyên được cung cấp trong repo này tương đối đầy đủ, theo nghĩa là chúng có thể được chạy cho các mẫu hiện có. Nếu bạn muốn giải các lưới của riêng mình, hãy lưu nó vào một tệp .json với nội dung có định dạng sau

    {
        "row": int,
        "column": int,
        "edges": list(list),
        "start": list,
        "target": list,
    }
Đây là từ điển chứa nội dung lưới của bạn, hãy xem thư mục 'Samples' để biết thêm thông tin. Nếu bạn không thích việc ghi và đọc các tệp trung gian, bạn có thể điều chỉnh code có sẵn để tạo và làm việc với dữ liệu bên trong chương trình.

### Lời cảm ơn
Nhóm xin chân thành cảm ơn các thầy và các bạn cùng lớp đã quan tâm và hỗ trợ trong suốt quá trình làm bài tập nhỏ này.