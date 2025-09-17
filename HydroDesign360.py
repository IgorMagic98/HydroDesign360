import tkinter as tk
from tkinter import messagebox, ttk
from sympy import symbols, Eq, solve

class DirectedGraphApp:
    NODE_RADIUS = 15

    def __init__(self, root):
        self.root = root
        self.root.title("Гидравлическая сеть с граничными условиями")

        style = ttk.Style()
        style.configure("Treeview", font=("Roboto Flex", 10))
        style.configure("Treeview.Heading", font=("Roboto Flex", 10, "bold"))

        # Канвас для визуализации графа
        self.canvas = tk.Canvas(root, width=700, height=500, bg='white')
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Привязка кликов мыши
        self.canvas.bind('<Button-1>', self.left_click_canvas)    # Левая кнопка
        self.canvas.bind('<Button-3>', self.right_click_canvas)   # Правая кнопка

        # Правая панель с вкладками
        right_frame = tk.Frame(root)
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        tab_control = ttk.Notebook(right_frame)
        tab_control.pack(fill=tk.BOTH, expand=True)

        # Вкладка с уравнениями и матрицей
        frame_eq = tk.Frame(tab_control)
        frame_eq.pack(fill=tk.BOTH, expand=True)

        self.matrix_text = tk.Text(frame_eq, width=40, height=10, font=("Roboto Flex", 10))
        self.matrix_text.pack(fill=tk.BOTH, expand=True)
        self.equations_text = tk.Text(frame_eq, width=40, height=10, font=("Roboto Flex", 10))
        self.equations_text.pack(fill=tk.BOTH, expand=True)
        self.solution_text = tk.Text(frame_eq, width=40, height=10, font=("Roboto Flex", 10))
        self.solution_text.pack(fill=tk.BOTH, expand=True)

        tab_control.add(frame_eq, text="Уравнения и матрица")

        # Вкладка с таблицей узлов (напор и тип)
        self.nodes_tree = ttk.Treeview(tab_control, columns=('H_val', 'y_val', 'Type'), show='headings')
        self.nodes_tree.heading('H_val', text='H (напор)')
        self.nodes_tree.column('H_val', width=100)
        self.nodes_tree.heading('y_val', text='y (высота)')
        self.nodes_tree.column('y_val', width=100)
        self.nodes_tree.heading('Type', text='Тип узла')
        self.nodes_tree.column('Type', width=100)
        tab_control.add(self.nodes_tree, text="Узлы")
        self.nodes_tree.bind("<Double-1>", self.on_node_edit)

        # Вкладка с таблицей рёбер
        self.edges_tree = ttk.Treeview(tab_control, columns=('Start', 'End', 'Q_val', 'K_val'), show='headings')
        self.edges_tree.heading('Start', text='Начало')
        self.edges_tree.column('Start', width=80)
        self.edges_tree.heading('End', text='Конец')
        self.edges_tree.column('End', width=80)
        self.edges_tree.heading('Q_val', text='Расход Q')
        self.edges_tree.column('Q_val', width=100)
        self.edges_tree.heading('K_val', text='Коэффициент K')
        self.edges_tree.column('K_val', width=130)
        tab_control.add(self.edges_tree, text="Рёбра")
        self.edges_tree.bind("<Double-1>", self.on_edge_edit)

        # Нижняя панель со списком входного и выходных узлов
        middle_frame = tk.Frame(right_frame, bg="lightgray", height=150)
        middle_frame.pack(fill=tk.BOTH, expand=False)

        input_frame = tk.Frame(middle_frame, bg="lightgray")
        input_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        tk.Label(input_frame, text="Входной узел", font=("Roboto Flex", 10)).pack()
        self.input_listbox = tk.Listbox(input_frame, font=("Roboto Flex", 10), height=1)
        self.input_listbox.pack(fill=tk.BOTH, expand=True)

        output_frame = tk.Frame(middle_frame, bg="lightgray")
        output_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        tk.Label(output_frame, text="Выходные узлы", font=("Roboto Flex", 10)).pack()
        self.output_listbox = tk.Listbox(output_frame, font=("Roboto Flex", 10))
        self.output_listbox.pack(fill=tk.BOTH, expand=True)

        # Нижняя панель кнопок
        btn_frame = tk.Frame(root)
        btn_frame.pack(side=tk.BOTTOM, fill=tk.X)
        tk.Button(btn_frame, text="Граничные условия", font=("Roboto Flex", 10), command=self.open_bc_dialog).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="Инициализация", font=("Roboto Flex", 10), command=self.find_paths_and_update_flows).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="Решить", font=("Roboto Flex", 10), command=self.recalculate).pack(side=tk.LEFT)

        # Основные данные графа
        self.nodes = []            # Список координат узлов [(x, y), ...]
        self.edges = []            # Список рёбер [(start_idx, end_idx), ...]
        self.H_values = {}         # Напоры H по узлам, ключ: индекс 1-based, значение: float
        self.node_heights = {}     # Высоты y по узлам, ключ: индекс 1-based, значение: float
        self.Q_values = {}         # Расходы по рёбрам, ключ: (start, end), значение: float
        self.K_values = {}         # Коэффициенты сопротивления K по рёбрам
        self.K = 1.0               # Значение K по умолчанию
        self.node_types = {}       # Типы узлов normal/input/output, ключ: индекс 0-based
        self.input_node = None
        self.output_nodes = set()

        self.all_paths = []
        self.equation_list = []
        self.selected_nodes = []

        self.redraw()
        self.update_trees()
        self.update_node_lists()

    def find_node(self, x, y):
        """Найти узел по координатам клика"""
        for idx, (nx, ny) in enumerate(self.nodes):
            dx = nx - x
            dy = ny - y
            dist_sq = dx * dx + dy * dy
            if dist_sq <= self.NODE_RADIUS * self.NODE_RADIUS:
                return idx
        return None

    def left_click_canvas(self, event):
        """Добавить узел или нарисовать ребро при последовательных кликах"""
        clicked = self.find_node(event.x, event.y)
        if clicked is None:
            self.add_node(event.x, event.y)
        else:
            self.selected_nodes.append(clicked)
            if len(self.selected_nodes) == 2:
                s, e = self.selected_nodes
                if s != e and (s, e) not in self.edges:
                    self.edges.append((s, e))
                self.selected_nodes.clear()
        self.redraw()
        self.update_trees()
        self.update_node_lists()

    def right_click_canvas(self, event):
        """Удалить узел или ребро по правому клику"""
        clicked_node = self.find_node(event.x, event.y)
        if clicked_node is not None:
            self.delete_node(clicked_node)
        else:
            for idx, (s, e) in enumerate(self.edges):
                if self.point_near_line(event.x, event.y, self.nodes[s][0], self.nodes[s][1], self.nodes[e][0], self.nodes[e][1]):
                    self.delete_edge(idx)
                    break
        self.redraw()
        self.update_trees()
        self.update_node_lists()

    def point_near_line(self, px, py, x1, y1, x2, y2, tol=5):
        """Проверяет, находится ли точка (px, py) на расстоянии tol от линии (x1,y1)-(x2,y2)"""
        from math import sqrt
        dx = x2 - x1
        dy = y2 - y1
        if dx == 0 and dy == 0:
            return sqrt((px - x1) ** 2 + (py - y1) ** 2) < tol
        t = max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / (dx * dx + dy * dy)))
        nearest_x = x1 + t * dx
        nearest_y = y1 + t * dy
        dist = sqrt((px - nearest_x) ** 2 + (py - nearest_y) ** 2)
        return dist < tol

    def add_node(self, x, y):
        """Добавить узел и пометить типом normal"""
        self.nodes.append((x, y))
        idx = len(self.nodes) - 1
        self.node_types[idx] = "normal"  # тип по умолчанию
        # Для нового узла установим напор и высоту по умолчанию
        # self.H_values[idx + 1] = 0.0
        self.node_heights[idx + 1] = 0.0
        h_val = self.H_values.get(idx + 1, 0.0)
        y_val = self.node_heights.get(idx + 1, 0.0)
        self.nodes_tree.insert('', 'end', iid=f"node{idx+1}",
                            values=(f"{h_val:.2f}", f"{y_val:.2f}", "normal"))


    def delete_node(self, node_idx):
        """Удалить узел и связанные с ним рёбра"""
        self.nodes.pop(node_idx)
        self.node_types.pop(node_idx, None)
        self.H_values = {k: v for k, v in self.H_values.items() if k != node_idx + 1}
        self.node_heights = {k: v for k, v in self.node_heights.items() if k != node_idx + 1}
        self.edges = [(s, e) for s, e in self.edges if s != node_idx and e != node_idx]
        self.output_nodes.discard(node_idx)
        if self.input_node == node_idx:
            self.input_node = None
        # Сдвиги индексов после удаления
        self.node_types = {k - 1 if k > node_idx else k: v for k, v in self.node_types.items()}
        if self.input_node is not None and self.input_node > node_idx:
            self.input_node -= 1
        self.output_nodes = {i - 1 if i > node_idx else i for i in self.output_nodes}
        self.update_trees()
        self.update_node_lists()
        self.redraw()
        self.nodes_tree.delete(f"node{node_idx+1}")

    def delete_edge(self, edge_idx):
        """Удалить ребро по индексу"""
        edge = self.edges.pop(edge_idx)
        self.Q_values.pop((edge[0] + 1, edge[1] + 1), None)
        self.update_trees()
        self.redraw()

    def open_bc_dialog(self):
        dlg = tk.Toplevel(self.root)
        dlg.title("Задать граничные условия")

        tk.Label(dlg, text="Напор (H) узла:").pack(padx=20)
        h_entry = tk.Entry(dlg)
        h_entry.pack(padx=20, pady=5)

        tk.Label(dlg, text="Высотная отметка (y) узла:").pack(padx=20)
        y_entry = tk.Entry(dlg)
        y_entry.pack(padx=20, pady=5)

        tk.Label(dlg, text="Тип граничных условий:").pack(pady=5)
        var_type = tk.StringVar(value="node")
        rb_node = tk.Radiobutton(dlg, text="Узел", variable=var_type, value="node")
        rb_edge = tk.Radiobutton(dlg, text="Ребро", variable=var_type, value="edge")
        rb_node.pack(anchor=tk.W, padx=20)
        rb_edge.pack(anchor=tk.W, padx=20)

        elements_var = tk.StringVar()
        elements_combobox = ttk.Combobox(dlg, state="readonly", textvariable=elements_var, width=30)
        elements_combobox.pack(pady=5, padx=20)

        type_node_var = tk.StringVar(value="normal")
        type_frame = tk.Frame(dlg)
        tk.Label(type_frame, text="Тип узла (если выбран узел):").pack(anchor=tk.W)
        tk.Radiobutton(type_frame, text="Входной", variable=type_node_var, value="input").pack(anchor=tk.W)
        tk.Radiobutton(type_frame, text="Выходной", variable=type_node_var, value="output").pack(anchor=tk.W)
        tk.Radiobutton(type_frame, text="Обычный узел", variable=type_node_var, value="normal").pack(anchor=tk.W)
        type_frame.pack(pady=5)

        k_label = None
        k_entry = None  # поле для ввода K у ребра

        def on_element_selected(event=None):
            t = var_type.get()
            if t != "node":
                return
            try:
                idx = int(elements_combobox.get().split()[1]) - 1
            except Exception:
                idx = 0

            h_val = self.H_values.get(idx + 1, 0.0)
            y_val = self.node_heights.get(idx + 1, 0.0)
            h_entry.config(state='normal')
            h_entry.delete(0, tk.END)
            h_entry.insert(0, str(h_val))
            y_entry.delete(0, tk.END)
            y_entry.insert(0, str(y_val))

            type_node_var.set(self.node_types.get(idx, 'normal'))

            if type_node_var.get() == "normal":
                h_entry.config(state='disabled')
            else:
                h_entry.config(state='normal')

        elements_combobox.bind("<<ComboboxSelected>>", on_element_selected)

        def on_type_node_changed(*args):
            if type_node_var.get() == "normal":
                h_entry.config(state='disabled')
            else:
                h_entry.config(state='normal')

        type_node_var.trace_add("write", on_type_node_changed)

        def update_elements(*args):
            nonlocal k_label, k_entry, h_entry, y_entry
            t = var_type.get()
            if t == "node":
                vals = [f"Узел {i+1}" for i in range(len(self.nodes))]
                elements_combobox['values'] = vals
                if vals:
                    current_selection = elements_combobox.get()
                    if current_selection not in vals:
                        elements_combobox.current(0)
                    else:
                        elements_combobox.set(current_selection)

                try:
                    idx = int(elements_combobox.get().split()[1]) - 1
                except Exception:
                    idx = 0
                h_val = self.H_values.get(idx + 1, 0.0)
                y_val = self.node_heights.get(idx + 1, 0.0)
                h_entry.config(state='normal')
                h_entry.delete(0, tk.END)
                h_entry.insert(0, str(h_val))
                y_entry.delete(0, tk.END)
                y_entry.insert(0, str(y_val))
                # НЕ вызываем type_node_var.set здесь, чтобы не сбрасывать выбор

                type_frame.pack(pady=5)

                if type_node_var.get() == "normal":
                    h_entry.config(state='disabled')
                else:
                    h_entry.config(state='normal')

                if k_entry:
                    k_entry.pack_forget()
                if k_label:
                    k_label.pack_forget()

                if not y_entry.winfo_ismapped():
                    y_entry.pack(padx=20, pady=5)
            else:
                vals = [f"Ребро {s+1} -> {e+1}" for s, e in self.edges]
                elements_combobox['values'] = vals
                if vals:
                    current_selection = elements_combobox.get()
                    if current_selection not in vals:
                        elements_combobox.current(0)
                    else:
                        elements_combobox.set(current_selection)

                type_frame.pack_forget()

                if h_entry.winfo_ismapped():
                    h_entry.pack_forget()
                if y_entry.winfo_ismapped():
                    y_entry.pack_forget()

                if not k_entry:
                    k_label = tk.Label(dlg, text="Коэффициент сопротивления K для ребра:")
                    k_label.pack(padx=20)
                    k_entry = tk.Entry(dlg)
                    k_entry.pack(padx=20, pady=5)
                else:
                    k_entry.pack(padx=20, pady=5)

        var_type.trace_add("write", update_elements)
        type_node_var.trace_add("write", update_elements)
        update_elements()

        def apply():
            t = var_type.get()
            el = elements_var.get()

            if t == "node":
                if not el:
                    messagebox.showerror("Ошибка", "Выберите узел")
                    return
                idx = int(el.split()[1]) - 1
                node_type = type_node_var.get()
                self.node_types[idx] = node_type

                try:
                    y_val = float(y_entry.get())
                    self.node_heights[idx + 1] = y_val
                except ValueError:
                    messagebox.showerror("Ошибка", "Высотная отметка y должна быть числом")
                    return

                if node_type == "normal":
                    # Напор H не задается, удаляем если был
                    self.H_values.pop(idx + 1, None)
                    if self.input_node == idx:
                        self.input_node = None
                    self.output_nodes.discard(idx)
                else:
                    try:
                        h_val = float(h_entry.get())
                        self.H_values[idx + 1] = h_val
                    except ValueError:
                        messagebox.showerror("Ошибка", "Напор H должен быть числом")
                        return

                    if node_type == "input":
                        self.input_node = idx
                        self.output_nodes.discard(idx)
                    elif node_type == "output":
                        try:
                            h_val = float(h_entry.get())
                            self.H_values[idx + 1] = h_val
                        except ValueError:
                            messagebox.showerror("Ошибка", "Напор H должен быть числом")
                            return
                        self.output_nodes.add(idx)
                        if self.input_node == idx:
                            self.input_node = None


            else:
                if not el:
                    messagebox.showerror("Ошибка", "Выберите ребро")
                    return
                parts = el.split()
                start = int(parts[1].split("->")[0]) - 1
                end = int(parts[2]) - 1
                val_str = val_entry.get()
                try:
                    val = float(val_str)
                except ValueError:
                    messagebox.showerror("Ошибка", "Значение должно быть числом")
                    return
                self.Q_values[(start + 1, end + 1)] = val
                try:
                    k_val = float(k_entry.get())
                except Exception:
                    k_val = self.K
                self.K_values[(start + 1, end + 1)] = k_val

            self.update_node_lists()
            self.update_trees()
            self.redraw()
            self.update_equations()
            self.recalculate()

        tk.Button(dlg, text="Применить", command=apply).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(dlg, text="Закрыть", command=dlg.destroy).pack(side=tk.LEFT, padx=5, pady=5)

    def update_node_lists(self):
        """Обновить списки входных и выходных узлов"""
        self.input_listbox.delete(0, tk.END)
        self.output_listbox.delete(0, tk.END)
        if self.input_node is not None:
            H = self.H_values.get(self.input_node + 1, 0.0)
            self.input_listbox.insert(tk.END, f"Вход: Узел {self.input_node + 1}, H={H:.2f}")
        for idx in sorted(self.output_nodes):
            H = self.H_values.get(idx + 1, 0.0)
            self.output_listbox.insert(tk.END, f"Выход: Узел {idx + 1}, H={H:.2f}")

    def find_all_paths(self, graph, start, end, path=None):
        """Рекурсивный поиск всех путей от start до end"""
        if path is None:
            path = []
        path = path + [start]
        if start == end:
            return [path]
        if start not in graph:
            return []
        paths = []
        for node in graph[start]:
            if node not in path:
                newpaths = self.find_all_paths(graph, node, end, path)
                for newpath in newpaths:
                    paths.append(newpath)
        return paths

    def find_paths_and_update_flows(self):
        """Инициализация и подсчет расходов по путям между входом и выходами"""
        if self.input_node is None or not self.output_nodes:
            messagebox.showwarning("Ошибка", "Задайте входной и выходные узлы")
            return

        graph = {i: [] for i in range(len(self.nodes))}
        for s, e in self.edges:
            graph[s].append(e)

        all_paths = []
        for out_node in self.output_nodes:
            all_paths.extend(self.find_all_paths(graph, self.input_node, out_node))

        if not all_paths:
            messagebox.showinfo("Инфо", "Пути не найдены")
            return

        self.all_paths = all_paths

        edges_in_paths = set()
        for path in all_paths:
            for i in range(len(path) - 1):
                edges_in_paths.add((path[i], path[i + 1]))

        # Наложение начальных значений расходов на рёбра
        self.Q_values = {k: 0.0 for k in self.Q_values.keys()}
        for s, e in self.edges:
            key = (s + 1, e + 1)
            if (s, e) in edges_in_paths:
                self.Q_values[key] = 1.0

        self.update_equations()
        self.recalculate()
        
        self.redraw()
        self.update_trees()
        self.update_node_lists()

    def update_equations(self):
        """Формировать систему уравнений по сохранению массы и балансу напоров"""

        self.equation_list.clear()
        self.equations_text.delete('1.0', tk.END)
        self.solution_text.delete('1.0', tk.END)
        self.equation_count = 0

        n = len(self.nodes)
        m = len(self.edges)

        # Матрица инцидентности (n узлов на m рёбер)
        matrix = [[0] * m for _ in range(n)]
        for col, (start, end) in enumerate(self.edges):
            matrix[start][col] = -1
            matrix[end][col] = 1

        # Исключаем узлы с <2 связей (будут зависимы)
        rows_to_keep = []
        for i in range(n):
            if sum(1 for val in matrix[i] if val != 0) >= 2:
                rows_to_keep.append(i)

        if not rows_to_keep:
            messagebox.showerror("Ошибка", "Нет узлов с двумя и более связями.")
            return

        reduced_matrix = [matrix[i] for i in rows_to_keep]

        Q_symbols = [symbols(f"Q{self.edges[j][0] + 1}{self.edges[j][1] + 1}") for j in range(m)]

        # Уравнения сохранения массы: сумма потоков в узле = 0
        for idx, row in enumerate(reduced_matrix):
            expr = sum(row[j] * Q_symbols[j] for j in range(m))
            q_syms = [sym for sym in expr.free_symbols if str(sym).startswith('Q')]
            if len(q_syms) < 1:
                continue
            self.equation_count += 1
            eq_name = f"eq{self.equation_count}"
            equation = Eq(expr, 0)
            self.equation_list.append(equation)

            terms = []
            for j, val in enumerate(row):
                if val == 0:
                    continue
                label = f"Q{self.edges[j][0] + 1}{self.edges[j][1] + 1}"
                sign = '+' if val > 0 else '-'
                terms.append(f"{sign}{label}")
            expr_str = "".join(terms)
            if expr_str.startswith('+'):
                expr_str = expr_str[1:]
            self.equations_text.insert(tk.END, f"{eq_name}: {expr_str} = 0\n")

        # Уравнения напорного баланса по рёбрам:
        for j, (s, e) in enumerate(self.edges):
            self.equation_count += 1
            eq_name = f"eq{self.equation_count}"
            Q = Q_symbols[j]

            # H, которые либо заданы (константы), либо символичны (переменные)
            Hs = self.H_values.get(s + 1, None)
            He = self.H_values.get(e + 1, None)
            Hs_sym = symbols(f"H{s + 1}")
            He_sym = symbols(f"H{e + 1}")
            H1 = Hs if Hs is not None else Hs_sym
            H2 = He if He is not None else He_sym

            K_val = self.K_values.get((s + 1, e + 1), self.K)
            Q_val = self.Q_values.get((s + 1, e + 1), 0.0)

            y_i = self.node_heights.get(s + 1, 0.0)  # высота узла начала ребра
            y_j = self.node_heights.get(e + 1, 0.0)  # высота узла конца ребра

            if abs(Q_val) < 1e-9:
                equation = Eq(Q, 0)
                self.equations_text.insert(tk.END, f"{eq_name}: Q{s + 1}{e + 1} = 0\n")
            else:
                # Уравнение напорного баланса с учетом разницы высот и сопротивления
                equation = Eq(H1 - H2 + (y_j - y_i) - K_val * Q * Q, 0)
                H1_str = f"{Hs:.2f}" if Hs is not None else f"H{s + 1}"
                H2_str = f"{He:.2f}" if He is not None else f"H{e + 1}"
                self.equations_text.insert(tk.END,
                    f"{eq_name}: {H1_str} - {H2_str} + ({y_j:.2f} - {y_i:.2f}) - {K_val:.2f}*Q{s + 1}{e + 1}^2 = 0\n")
            self.equation_list.append(equation)

        # Подсчет неизвестных
        all_vars = set()
        for eq in self.equation_list:
            all_vars.update(eq.free_symbols)
        known_H = {symbols(f"H{i}") for i in self.H_values.keys()}
        unknown_vars = {v for v in all_vars if v not in known_H}

        self.equations_text.insert(tk.END, f"\nЧисло уравнений: {len(self.equation_list)}\n")
        self.equations_text.insert(tk.END, f"Число неизвестных: {len(unknown_vars)}\n")

        # Вывести граничные условия - напоры и высоты
        bc_text = "Граничные условия:\n"
        for i in range(len(self.nodes)):
            h_val = self.H_values.get(i + 1, None)
            y_val = self.node_heights.get(i + 1, 0.0)
            # if h_val is not None:
            if h_val is not None and h_val != 0.0:
                bc_text += f"Узел {i+1}: H={h_val:.2f}, y={y_val:.2f}\n"
            else:
                bc_text += f"Узел {i+1}: H=не задан, y={y_val:.2f}\n"
        self.solution_text.insert(tk.END, bc_text)

        # Вывести расходы по рёбрам
        self.solution_text.insert(tk.END, "\nРасходы на рёбрах:\n")
        for s, e in self.edges:
            key = (s + 1, e + 1)
            q_val = self.Q_values.get(key, 0.0)
            self.solution_text.insert(tk.END, f"Ребро {key}: Q={q_val:.2f}\n")

    def recalculate(self):
        """Решение системы уравнений с подстановкой известных значений"""
        self.equations_text.delete('1.0', tk.END)
        self.solution_text.delete('1.0', tk.END)
        self.update_equations()
        try:
            variables = []
            for eq in self.equation_list:
                for symb in eq.free_symbols:
                    if symb not in variables:
                        variables.append(symb)

            # Подстановка известных напоров H для решения
            subs = {}
            for idx, val in self.H_values.items():
                sym = symbols(f"H{idx}")
                subs[sym] = val

            equations_to_solve = [eq.subs(subs) for eq in self.equation_list]
            vars_to_solve = [v for v in variables if v not in subs]

            solution = solve(equations_to_solve, vars_to_solve, dict=True)

            if solution:
                sol = solution[0]
                self.last_solution = sol
                for var, val in sol.items():
                    name = str(var)
                    if name.startswith("Q"):
                        key = (int(name[1]), int(name[2]))
                        self.Q_values[key] = float(val.evalf())
                    elif name.startswith("H"):
                        key = int(name[1:])
                        self.H_values[key] = float(val.evalf())
                self.solution_text.insert(tk.END, "\nРешение системы уравнений:\n")
                for var, val in sol.items():
                    self.solution_text.insert(tk.END, f"{var} = {val.evalf():.4f}\n")
            else:
                self.solution_text.insert(tk.END, "\nРешение не найдено или система несовместима.\n")

        except Exception as e:
            self.solution_text.insert(tk.END, f"\nОшибка при решении: {e}\n")

        self.redraw()
        self.update_node_lists()
        self.update_trees()
        self.update_equations()

    def redraw(self):
        """Перерисовка графа на canvas"""
        self.canvas.delete("all")

        edges_in_paths = set()
        for p in self.all_paths:
            for i in range(len(p) - 1):
                edges_in_paths.add((p[i], p[i + 1]))

        # Рисуем рёбра с цветом, указывающим на принадлежность к найденным путям
        for idx, (s, e) in enumerate(self.edges, start=1):
            x1, y1 = self.nodes[s]
            x2, y2 = self.nodes[e]
            color = "black"
            width = 2
            if (s, e) in edges_in_paths:
                color = "green"
                width = 3
            self.draw_arrow(x1, y1, x2, y2, color=color, width=width)
            q_val = self.Q_values.get((s + 1, e + 1), None)
            if q_val is not None:
                self.canvas.create_text((x1 + x2) / 2, (y1 + y2) / 2 + 10,
                                        text=f"Q={q_val:.2f}", fill='blue', font=("Roboto Flex", 10))

        # Рисуем узлы с цветом по типу
        for idx, (x, y) in enumerate(self.nodes, start=1):
            node_type = self.node_types.get(idx - 1, "normal")
            fill_color = "lightblue"
            if node_type == "input":
                fill_color = "lightgreen"
            elif node_type == "output":
                fill_color = "lightcoral"
            self.canvas.create_oval(x - self.NODE_RADIUS, y - self.NODE_RADIUS,
                                    x + self.NODE_RADIUS, y + self.NODE_RADIUS,
                                    fill=fill_color, outline="black")
            self.canvas.create_text(x, y, text=str(idx), font=("Roboto Flex", 12, "bold"))

            # Отображение высоты y под узлом
            y_val = self.node_heights.get(idx, 0.0)
            self.canvas.create_text(x, y + self.NODE_RADIUS + 10,
                                    text=f"y={y_val:.2f}", font=("Roboto Flex", 10), fill="black")

            # Отображение напора H над узлом (если задан)
            h_val = self.H_values.get(idx, None)
            if h_val is not None:
                self.canvas.create_text(x, y - self.NODE_RADIUS - 10,
                                        text=f"H={h_val:.2f}", font=("Roboto Flex", 10), fill="black")

        self.update_incidence_matrix()
        self.update_node_lists()
        self.update_trees()

    def draw_arrow(self, x1, y1, x2, y2, color="black", width=2):
        """Рисует стрелку с учетом радиуса узла"""
        from math import atan2, cos, sin
        dx = x2 - x1
        dy = y2 - y1
        angle = atan2(dy, dx)
        r = self.NODE_RADIUS
        start_x = x1 + r * cos(angle)
        start_y = y1 + r * sin(angle)
        end_x = x2 - r * cos(angle)
        end_y = y2 - r * sin(angle)
        self.canvas.create_line(start_x, start_y, end_x, end_y, arrow=tk.LAST, width=width, fill=color)

    def update_incidence_matrix(self):
        """Обновить матрицу инцидентности и отобразить ее в текстовом поле"""
        n = len(self.nodes)
        real_edges = list(self.edges)  # Текущие ребра
        m = len(real_edges)

        matrix = [[0] * m for _ in range(n)]
        for col, (start, end) in enumerate(real_edges):
            matrix[start][col] = -1
            matrix[end][col] = 1

        rows_to_keep = []
        for i in range(n):
            if sum(1 for v in matrix[i] if v != 0) >= 2:
                rows_to_keep.append(i)

        text = "Матрица инцидентности (узлы x рёбра):\n\n"
        header = ["e" + str(i + 1) for i in range(m)]
        text += "\t" + "\t".join(header) + "\n"
        for i in rows_to_keep:
            text += "n" + str(i + 1) + "\t" + "\t".join(map(str, matrix[i])) + "\n"

        self.matrix_text.delete(1.0, tk.END)
        self.matrix_text.insert(tk.END, text)

    def update_trees(self):
        """Обновить таблицы узлов и рёбер"""
        self.nodes_tree.delete(*self.nodes_tree.get_children())
        for idx, (x, y) in enumerate(self.nodes, start=1):
            h_val = self.H_values.get(idx, "")
            y_val = self.node_heights.get(idx, "")
            node_type = self.node_types.get(idx - 1, "normal")
            if isinstance(h_val, (int, float)):
                h_val = f"{h_val:.2f}"
            if isinstance(y_val, (int, float)):
                y_val = f"{y_val:.2f}"
            self.nodes_tree.insert('', 'end', iid=f"node{idx}", values=(h_val, y_val, node_type))

        self.edges_tree.delete(*self.edges_tree.get_children())
        for idx, (s, e) in enumerate(self.edges, start=1):
            q_val = self.Q_values.get((s + 1, e + 1), 0.0)
            k_val = self.K_values.get((s + 1, e + 1), self.K)
            self.edges_tree.insert('', 'end', iid=f"edge{idx}",
                                   values=(s + 1, e + 1, f"{q_val:.2f}", f"{k_val:.2f}"))

    def on_node_edit(self, event):
        """Обработка двойного щелчка в таблице узлов"""
        self.edit_cell(event, self.nodes_tree, self.update_node_value)

    def on_edge_edit(self, event):
        """Обработка двойного щелчка в таблице рёбер"""
        self.edit_cell(event, self.edges_tree, self.update_edge_value)

    def edit_cell(self, event, tree, callback):
        """Реализация редактирования ячейки в таблицах"""
        x, y, widget = event.x, event.y, event.widget
        rowid = widget.identify_row(y)
        column = widget.identify_column(x)
        if not rowid or not column:
            return
        col_num = int(column[1:]) - 1
        bbox = widget.bbox(rowid, column)
        if not bbox:
            return
        xabs, yabs, width, height = bbox
        entry = tk.Entry(widget)
        entry.place(x=xabs, y=yabs, width=width, height=height)
        value = widget.set(rowid, column)
        entry.insert(0, value)
        entry.focus_set()

        def on_enter(event):
            new_val = entry.get()
            try:
                new_val_f = float(new_val)
            except ValueError:
                messagebox.showerror("Ошибка", "Введите число")
                return
            callback(rowid, col_num, new_val_f)
            entry.destroy()

        def on_focus_out(event):
            entry.destroy()

        entry.bind("<Return>", on_enter)
        entry.bind("<FocusOut>", on_focus_out)

    def update_node_value(self, rowid, col_num, value):
        """Обновление данных узла после редактирования в таблице"""
        node_idx = int(rowid.replace("node", "")) - 1
        if col_num == 0:
            # Обновляем напор H
            self.H_values[node_idx + 1] = value
        elif col_num == 1:
            # Обновляем высоту y
            self.node_heights[node_idx + 1] = value
        elif col_num == 2:
            # Обновляем тип узла ('normal', 'input', 'output')
            str_value = str(value)
            if str_value in ("normal", "input", "output"):
                self.node_types[node_idx] = str_value
                if str_value == "input":
                    self.input_node = node_idx
                    self.output_nodes.discard(node_idx)
                elif str_value == "output":
                    self.output_nodes.add(node_idx)
                    if self.input_node == node_idx:
                        self.input_node = None
                else:
                    if self.input_node == node_idx:
                        self.input_node = None
                    self.output_nodes.discard(node_idx)
        self.update_trees()
        self.recalculate()

    def update_edge_value(self, rowid, col_num, value):
        edge_idx = int(rowid.replace("edge", "")) - 1
        edge = self.edges[edge_idx]
        if col_num == 2:  # Изменение расхода Q
            self.Q_values[(edge[0] + 1, edge[1] + 1)] = value
        elif col_num == 3:  # Изменение коэффициента K
            self.K_values[(edge[0] + 1, edge[1] + 1)] = value
            # Пересчитать цепочки и уравнения после изменения K
            self.find_paths_and_update_flows()
            self.recalculate()


if __name__ == "__main__":
    root = tk.Tk()
    app = DirectedGraphApp(root)
    root.mainloop() #ghfdgh
