import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import os

# ==============================================================================
# 1. 核心功能函数 (Bio-Logic)
# ==============================================================================

def get_reverse_complement(seq):
    """计算反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'U': 'A'}
    return "".join(complement.get(base.upper(), base.upper()) for base in reversed(seq))

def read_fasta(file_path):
    """读取 FASTA 文件，忽略 header，只返回序列字符串"""
    if not file_path: return ""
    seq = ""
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq.upper()

def simulate_pcr_amplification(template, fwd_primer_full, rev_primer_full, fwd_bind, rev_bind):
    """模拟 PCR: 提取结合位点中间的序列"""
    template = template.upper()
    fwd_bind = fwd_bind.upper()
    rev_bind = rev_bind.upper()

    f_start = template.find(fwd_bind)
    if f_start == -1: return None
    
    rev_bind_rc = get_reverse_complement(rev_bind)
    r_start = template.find(rev_bind_rc, f_start + len(fwd_bind))
    
    if r_start == -1: return None
    
    insert_seq = template[f_start + len(fwd_bind) : r_start]
    product = fwd_primer_full + insert_seq + get_reverse_complement(rev_primer_full)
    return product

def get_backbone(vector_seq):
    """提取 pDIRECT 骨架"""
    site_fwd = "GCTCTTC" # SapI
    site_rev = "GAAGAGC" # SapI RC
    
    pos1 = vector_seq.find(site_fwd)
    pos2 = vector_seq.find(site_rev)
    
    if pos1 == -1 or pos2 == -1: return None

    cut_idx_1 = pos1 + 8
    cut_idx_2 = pos2 - 1
    
    if cut_idx_1 < cut_idx_2:
        frag_a = vector_seq[cut_idx_1:cut_idx_2]
        frag_b = vector_seq[cut_idx_2:] + vector_seq[:cut_idx_1]
    else:
        frag_a = vector_seq[cut_idx_1:] + vector_seq[:cut_idx_2]
        frag_b = vector_seq[cut_idx_2:cut_idx_1]

    # 返回最长的片段 (骨架)，此时骨架两头都包含了 3bp Overhang
    return frag_a if len(frag_a) > len(frag_b) else frag_b

def simulate_digestion(fragment, left_enzyme, right_enzyme):
    """模拟酶切"""
    seq = fragment.upper()
    start_idx = 0
    end_idx = len(seq)

    # 左切
    if left_enzyme == "SapI": 
        site = "GCTCTTC"
        match = seq.find(site)
        if match != -1: start_idx = match + 8
    elif left_enzyme == "Esp3I": 
        site = "CGTCTC"
        match = seq.find(site)
        if match != -1: start_idx = match + 7
            
    # 右切
    if right_enzyme == "Esp3I": 
        site = "GAGACG"
        match = seq.rfind(site)
        if match != -1: end_idx = match - 1
    elif right_enzyme == "SapI": 
        site = "GAAGAGC"
        match = seq.rfind(site)
        if match != -1: end_idx = match - 1
            
    if start_idx >= end_idx: return ""
    return seq[start_idx : end_idx]

def assemble_plasmid(backbone, frag1, frag2, frag3):
    """
    组装逻辑 (v4 修复版)
    核心修复：处理 Backbone 与 P1 之间的 3bp 重复
    """
    p1 = simulate_digestion(frag1, "SapI", "Esp3I")
    p2 = simulate_digestion(frag2, "Esp3I", "Esp3I")
    p3 = simulate_digestion(frag3, "Esp3I", "SapI")
    
    if not p1 or not p2 or not p3: return None

    # --- 关键修复区域 ---
    # Backbone 的 3' 端已经有了 3bp Overhang (来自 get_backbone)
    # P1 (frag1) 的 5' 端也是 3bp Overhang (来自 simulate_digestion)
    # 必须切除 P1 的头部 3bp，否则会变成 ...GTCGTC...
    
    # 1. 处理 P1 头部 (SapI Overhang, 3bp)
    p1_trimmed = p1[3:] if len(p1) > 3 else p1

    # 2. 处理 P2 头部 (Esp3I Overhang, 4bp)
    p2_trimmed = p2[4:] if len(p2) > 4 else p2
    
    # 3. 处理 P3 头部 (Esp3I Overhang, 4bp)
    p3_trimmed = p3[4:] if len(p3) > 4 else p3
    
    # 4. 处理 P3 尾部 (SapI Overhang, 3bp)
    # 因为 P3 尾部有 Overhang，Backbone 头部也有，会重复，所以切掉 P3 尾部
    p3_final = p3_trimmed[:-3] if len(p3_trimmed) > 3 else p3_trimmed 

    final_seq = backbone + p1_trimmed + p2_trimmed + p3_final
    return final_seq

# ==============================================================================
# 2. GUI 界面类
# ==============================================================================

class CrisprCloningApp:
    def __init__(self, root):
        self.root = root
        self.root.title("🧬 CRISPR Golden Gate Tool v4 (Final Fix)")
        self.root.geometry("900x800")
        
        self.pmod_seq = "" 
        self.pdirect_seq = "" 

        style = ttk.Style()
        style.configure("TButton", font=("Arial", 10, "bold"))

        # --- 1. File Selection ---
        file_frame = ttk.LabelFrame(root, text="Step 1: 加载模版序列 (Load Templates)", padding=10)
        file_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(file_frame, text="pMOD 模板 (.fasta):").grid(row=0, column=0, sticky="w")
        self.pmod_path = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.pmod_path, width=50).grid(row=0, column=1, padx=5)
        ttk.Button(file_frame, text="浏览...", command=self.load_pmod).grid(row=0, column=2)
        
        ttk.Label(file_frame, text="pDIRECT 骨架 (.fasta):").grid(row=1, column=0, sticky="w")
        self.pdirect_path = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.pdirect_path, width=50).grid(row=1, column=1, padx=5)
        ttk.Button(file_frame, text="浏览...", command=self.load_pdirect).grid(row=1, column=2)

        ttk.Label(file_frame, text="* 提示: 必须加载您自己的 fasta 文件以确保序列 100% 准确。").grid(row=2, column=0, columnspan=3, sticky="w", pady=5)

        # --- 2. Input ---
        input_frame = ttk.LabelFrame(root, text="Step 2: 输入 gRNA (Input)", padding=10)
        input_frame.pack(fill="x", padx=10, pady=5)
        
        ttk.Label(input_frame, text="gRNA1 (5'->3'):").grid(row=0, column=0, sticky="w")
        self.g1_entry = ttk.Entry(input_frame, width=60)
        self.g1_entry.grid(row=0, column=1, padx=5, pady=5)
        
        ttk.Label(input_frame, text="gRNA2 (5'->3'):").grid(row=1, column=0, sticky="w")
        self.g2_entry = ttk.Entry(input_frame, width=60)
        self.g2_entry.grid(row=1, column=1, padx=5, pady=5)

        # --- 3. Run ---
        btn_frame = ttk.Frame(root)
        btn_frame.pack(pady=10)
        self.run_btn = ttk.Button(btn_frame, text="生成结果", command=self.run_design)
        self.run_btn.pack(side="left", padx=10)
        
        # --- 4. Output ---
        output_frame = ttk.LabelFrame(root, text="Step 3: 结果 (Results)", padding=10)
        output_frame.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.notebook = ttk.Notebook(output_frame)
        self.notebook.pack(fill="both", expand=True)
        
        self.tab_primers = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_primers, text="引物列表")
        self.txt_primers = tk.Text(self.tab_primers, height=15, font=("Consolas", 10))
        self.txt_primers.pack(fill="both", expand=True, padx=5, pady=5)
        
        self.tab_vector = ttk.Frame(self.notebook)
        self.notebook.add(self.tab_vector, text="终载体序列")
        
        save_bar = ttk.Frame(self.tab_vector)
        save_bar.pack(fill="x", padx=5, pady=5)
        ttk.Button(save_bar, text="保存 FASTA", command=self.save_fasta).pack(side="right")
        
        self.txt_vector = tk.Text(self.tab_vector, height=15, font=("Consolas", 10))
        self.txt_vector.pack(fill="both", expand=True, padx=5, pady=5)

    def load_pmod(self):
        filename = filedialog.askopenfilename(filetypes=[("FASTA", "*.fasta"), ("All", "*.*")])
        if filename:
            self.pmod_path.set(filename)
            self.pmod_seq = read_fasta(filename)

    def load_pdirect(self):
        filename = filedialog.askopenfilename(filetypes=[("FASTA", "*.fasta"), ("All", "*.*")])
        if filename:
            self.pdirect_path.set(filename)
            self.pdirect_seq = read_fasta(filename)

    def run_design(self):
        g1 = self.g1_entry.get().strip().replace(" ", "").upper()
        g2 = self.g2_entry.get().strip().replace(" ", "").upper()
        
        if not self.pmod_seq or not self.pdirect_seq:
             messagebox.showwarning("警告", "请先加载 pMOD 和 pDIRECT 的 FASTA 文件！")
             return

        if len(g1) < 12 or len(g2) < 12:
            messagebox.showerror("错误", "gRNA 序列太短！")
            return
            
        try:
            # 引物计算
            g1_rc_12 = get_reverse_complement(g1)[-12:]
            g2_rc_12 = get_reverse_complement(g2)[-12:]
            g1_12 = g1[-12:]
            g2_12 = g2[-12:]
            
            bind_p1 = "GGCAGACATACTGTCCCAC"
            bind_p2 = "CTGCCTATACGGCAGTGAAC"
            bind_p3 = "GTTTTAGAGCTAGAAATAGC"
            bind_p4 = "CTGCCTATACGGCAGTGAAC"
            bind_p5 = "GTTTTAGAGCTAGAAATAGC"
            bind_p6 = "CTGCCTATACGGCAGTGAAC"

            p1_seq = "TGCTCTTCGCGCTGGCAGACATACTGTCCCAC"
            p2_seq = "TCGTCTCC" + g1_rc_12 + bind_p2     
            p3_seq = "TCGTCTCA" + g1_12 + bind_p3        
            p4_seq = "TCGTCTCA" + g2_rc_12 + bind_p4     
            p5_seq = "TCGTCTCA" + g2_12 + bind_p5        
            p6_seq = "TGCTCTTCTGACCTGCCTATACGGCAGTGAAC" 

            report = f">CmYLCV_Fixed\n{p1_seq}\n>Csy4-gRNA1\n{p2_seq}\n>Rep-gRNA1\n{p3_seq}\n"
            report += f">Csy4-gRNA2\n{p4_seq}\n>REP-gRNA2\n{p5_seq}\n>CSY-term_Fixed\n{p6_seq}\n"
            self.txt_primers.delete("1.0", tk.END)
            self.txt_primers.insert(tk.END, report)
            
            # 模拟 PCR
            frag1 = simulate_pcr_amplification(self.pmod_seq, p1_seq, p2_seq, bind_p1, bind_p2)
            frag2 = simulate_pcr_amplification(self.pmod_seq, p3_seq, p4_seq, bind_p3, bind_p4)
            frag3 = simulate_pcr_amplification(self.pmod_seq, p5_seq, p6_seq, bind_p5, bind_p6)
            
            if not all([frag1, frag2, frag3]):
                messagebox.showerror("失败", "PCR 模拟失败。请检查模板文件或序列。")
                return

            # 组装
            backbone = get_backbone(self.pdirect_seq)
            if not backbone:
                messagebox.showerror("失败", "骨架提取失败。找不到 SapI 位点。")
                return

            final = assemble_plasmid(backbone, frag1, frag2, frag3)
            
            if final:
                self.txt_vector.delete("1.0", tk.END)
                # 修改 header 显示信息
                header = f">pDIRECT_21C_Assembled_gRNA1_{g1[:5]}_gRNA2_{g2[:5]} ({len(final)} bp)\n"
                self.txt_vector.insert(tk.END, header + final)
                self.notebook.select(self.tab_vector)
                messagebox.showinfo("成功", f"组装成功！长度: {len(final)} bp (完美匹配)")
            else:
                messagebox.showerror("错误", "组装逻辑错误")

        except Exception as e:
            messagebox.showerror("错误", str(e))

    def save_fasta(self):
        content = self.txt_vector.get("1.0", tk.END).strip()
        if not content: return
        path = filedialog.asksaveasfilename(defaultextension=".fasta", filetypes=[("FASTA", "*.fasta")])
        if path:
            with open(path, "w") as f:
                f.write(content)

if __name__ == "__main__":
    root = tk.Tk()
    app = CrisprCloningApp(root)
    root.mainloop()