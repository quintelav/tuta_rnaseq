#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SAQE - Gerador de Tabelas Académicas e Científicas de Alta Resolução
Lê a planilha anotada da metanálise e exporta tabelas formatadas em LaTeX (booktabs),
Markdown de relatórios e um Dashboard HTML interativo responsivo.
"""

import os
import json
import re
import pandas as pd
import numpy as np

def escape_latex(text):
    """
    Escapa de forma estrita caracteres especiais do LaTeX para garantir
    a compilação bem-sucedida do documento sem quebras de código.
    """
    if not isinstance(text, str):
        return ""
    # Protege formatações científicas comuns e formata termos biológicos em itálico
    text = text.replace("Tuta absoluta", r"\textit{Tuta absoluta}")
    text = text.replace("Phthorimaea absoluta", r"\textit{Phthorimaea absoluta}")
    
    replacements = {
        '&': r'\&',
        '%': r'\%',
        '$': r'\$',
        '#': r'\#',
        '_': r'\_',
        '{': r'\{',
        '}': r'\}',
        '~': r'\textasciitilde{}',
        '^': r'\textasciicircum{}',
        '\\': r'\textbackslash{}',
    }
    for key, val in replacements.items():
        if key in text:
            # Garante que não duplica se já estiver escapado
            text = re.sub(r'(?<!\\)' + re.escape(key), val, text)
    return text

def gerar_tabela_latex(df_dados, caminho_saida):
    """
    Gera uma tabela LaTeX usando o pacote professional 'booktabs' (sem linhas verticais).
    Ideal para submissão direta em periódicos internacionais.
    """
    print("✍️  A gerar tabela profissional em LaTeX (estilo booktabs)...")
    
    linhas = []
    linhas.append(r"\begin{table}[htbp]")
    linhas.append(r"  \centering")
    linhas.append(r"  \caption{Mecanismo central de detoxificação de \textit{Phthorimaea absoluta} sob stresse de inseticidas: Principais transcritos diferencialmente ativos identificados via On-Score e On-Ratio (Metanálise Bono, 2021).}")
    linhas.append(r"  \label{tab:metanalise_detox_candidatos}")
    linhas.append(r"  \small")
    
    # Definição das colunas: ID, Anotação, Log2FC Abamectina, Spinosad, Emamectina, On-Score (UP), On-Ratio
    linhas.append(r"  \begin{tabular}{l p{6cm} c c c c c}")
    linhas.append(r"    \toprule")
    linhas.append(r"    \textbf{Transcript ID} & \textbf{Nome/Anotação Funcional} & \textbf{Abam.} & \textbf{Spin.} & \textbf{Emam.} & \textbf{On-Score} & \textbf{On-Ratio} \\")
    linhas.append(r"    & & \textbf{($\log_2$ FC)} & \textbf{($\log_2$ FC)} & \textbf{($\log_2$ FC)} & \textbf{(UP)} & \textbf{(UP)} \\")
    linhas.append(r"    \midrule")
    
    for idx, row in df_dados.iterrows():
        id_esc = escape_latex(str(idx))
        
        # Limpa anotações muito compridas para caberem de forma limpa na célula da tabela
        anot = str(row.get('Annotation', 'Sem descrição'))
        if len(anot) > 65:
            anot = anot[:62] + "..."
        anot_esc = escape_latex(anot)
        
        fc_abam = f"{row.get('Log2Ratio_abamectina', 0.0):.2f}" if 'Log2Ratio_abamectina' in df_dados.columns else "0.00"
        fc_spin = f"{row.get('Log2Ratio_spinosad', 0.0):.2f}" if 'Log2Ratio_spinosad' in df_dados.columns else "0.00"
        fc_emam = f"{row.get('Log2Ratio_emamectina', 0.0):.2f}" if 'Log2Ratio_emamectina' in df_dados.columns else "0.00"
        
        score = int(row.get('Bono_UP_Score', 0))
        ratio = f"{row.get('On_Ratio_UP', 0.0):.2f}"
        
        # Destaca em negrito na tabela os genes com On-Score igual a 3 (resposta universal)
        if score == 3:
            linhas.append(f"    \\textbf{{{id_esc}}} & \\textbf{{{anot_esc}}} & \\textbf{{{fc_abam}}} & \\textbf{{{fc_spin}}} & \\textbf{{{fc_emam}}} & \\textbf{{{score}}} & \\textbf{{{ratio}}} \\\\")
        else:
            linhas.append(f"    {id_esc} & {anot_esc} & {fc_abam} & {fc_spin} & {fc_emam} & {score} & {ratio} \\\\")
            
    linhas.append(r"    \bottomrule")
    linhas.append(r"  \end{tabular}")
    linhas.append(r"\end{table}")
    
    with open(caminho_saida, 'w', encoding='utf-8') as f:
        f.write("\n".join(linhas))
        
    print(f"   ✅ LaTeX guardado em: '{caminho_saida}'")

def gerar_tabela_markdown(df_dados, caminho_saida):
    """
    Gera uma tabela legível em formato Markdown para relatórios rápidos.
    """
    print("✍️  A gerar tabela de visualização rápida em Markdown...")
    
    linhas = []
    linhas.append("# Metanálise de Transcritoma - Tabela de Candidatos de Detoxificação (On-Score >= 2)")
    linhas.append("")
    linhas.append("| Transcript ID | Anotação Funcional / Gene | Família de Detox | Log2FC Abamectina | Log2FC Spinosad | Log2FC Emamectina | On-Score (UP) | On-Ratio |")
    linhas.append("|---|---|---|---|---|---|---|---|")
    
    for idx, row in df_dados.iterrows():
        anot = str(row.get('Annotation', 'Sem descrição'))
        if len(anot) > 65:
            anot = anot[:62] + "..."
            
        fam = str(row.get('Detox_Family', 'Outros'))
        
        fc_abam = f"{row.get('Log2Ratio_abamectina', 0.0):.2f}" if 'Log2Ratio_abamectina' in df_dados.columns else "0.00"
        fc_spin = f"{row.get('Log2Ratio_spinosad', 0.0):.2f}" if 'Log2Ratio_spinosad' in df_dados.columns else "0.00"
        fc_emam = f"{row.get('Log2Ratio_emamectina', 0.0):.2f}" if 'Log2Ratio_emamectina' in df_dados.columns else "0.00"
        
        score = int(row.get('Bono_UP_Score', 0))
        ratio = f"{row.get('On_Ratio_UP', 0.0):.2f}"
        
        linhas.append(f"| `{idx}` | {anot} | {fam} | {fc_abam} | {fc_spin} | {fc_emam} | **{score}** | {ratio} |")
        
    with open(caminho_saida, 'w', encoding='utf-8') as f:
        f.write("\n".join(linhas))
        
    print(f"   ✅ Markdown guardado em: '{caminho_saida}'")

def gerar_dashboard_html(df_dados, caminho_saida):
    """
    Gera um Dashboard HTML interativo completo, estilizado com Tailwind CSS.
    Injeta os dados em JSON para garantir funcionamento offline em qualquer navegador.
    """
    print("✍️  A compilar o Dashboard Interativo em HTML...")
    
    # Converte o DataFrame filtrado para uma estrutura JSON segura
    registros = []
    for idx, row in df_dados.iterrows():
        registros.append({
            'id': str(idx),
            'annotation': str(row.get('Annotation', 'Sem descrição molecular')),
            'family': str(row.get('Detox_Family', 'Outros')),
            'fc_abam': round(float(row.get('Log2Ratio_abamectina', 0.0)), 2) if 'Log2Ratio_abamectina' in df_dados.columns else 0.0,
            'fc_spin': round(float(row.get('Log2Ratio_spinosad', 0.0)), 2) if 'Log2Ratio_spinosad' in df_dados.columns else 0.0,
            'fc_emam': round(float(row.get('Log2Ratio_emamectina', 0.0)), 2) if 'Log2Ratio_emamectina' in df_dados.columns else 0.0,
            'score': int(row.get('Bono_UP_Score', 0)),
            'ratio': round(float(row.get('On_Ratio_UP', 0.0)), 2)
        })
        
    dados_json = json.dumps(registros, indent=2)
    
    html_content = f"""<!DOCTYPE html>
<html lang="pt-PT">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SAQE - Explorador Interativo de Metanálise (Bono, 2021)</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    <style>
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
        body {{
            font-family: 'Inter', sans-serif;
        }}
    </style>
</head>
<body class="bg-slate-50 text-slate-800 min-h-screen flex flex-col">

    <!-- CABEÇALHO -->
    <header class="bg-gradient-to-r from-emerald-700 to-teal-800 text-white shadow-md">
        <div class="max-w-7xl mx-auto px-4 py-6 sm:px-6 lg:px-8 flex flex-col sm:flex-row justify-between items-center">
            <div>
                <h1 class="text-2xl font-bold tracking-tight">SAQE - Explorador Metanalítico</h1>
                <p class="text-emerald-100 text-sm mt-1">Metanálise Transcriptómica de <span class="italic font-semibold">Phthorimaea absoluta</span> (Bono, 2021)</p>
            </div>
            <div class="mt-4 sm:mt-0 flex gap-3 text-xs bg-emerald-900/40 p-2 rounded-lg border border-emerald-500/20">
                <div class="text-right">
                    <span class="block text-emerald-300 font-medium">Desenho de Estudo</span>
                    <span class="text-white">3 Inseticidas (Spinosad, Abamectina, Emamectina)</span>
                </div>
            </div>
        </div>
    </header>

    <!-- PAINEL PRINCIPAL -->
    <main class="flex-grow max-w-7xl w-full mx-auto px-4 py-8 sm:px-6 lg:px-8">
        
        <!-- CARD METRICS -->
        <div class="grid grid-cols-1 md:grid-cols-4 gap-4 mb-8">
            <div class="bg-white p-5 rounded-xl shadow-sm border border-slate-200/80">
                <span class="text-xs text-slate-400 font-semibold uppercase tracking-wider block">Transcritos Candidatos</span>
                <span class="text-3xl font-bold text-slate-800 block mt-1" id="total-count">0</span>
                <span class="text-xs text-emerald-600 font-medium mt-1 inline-block"><i class="fa-solid fa-filter"></i> On-Score &ge; 2</span>
            </div>
            <div class="bg-white p-5 rounded-xl shadow-sm border border-slate-200/80">
                <span class="text-xs text-slate-400 font-semibold uppercase tracking-wider block">Frequência Core (Score = 3)</span>
                <span class="text-3xl font-bold text-indigo-600 block mt-1" id="core-count">0</span>
                <span class="text-xs text-slate-500 block mt-1">Ativo nos 3 tratamentos</span>
            </div>
            <div class="bg-white p-5 rounded-xl shadow-sm border border-slate-200/80">
                <span class="text-xs text-slate-400 font-semibold uppercase tracking-wider block">Alvos de Fase I (CYP/CCE)</span>
                <span class="text-3xl font-bold text-emerald-600 block mt-1" id="phase1-count">0</span>
                <span class="text-xs text-slate-500 block mt-1">Oxidação e hidrólise</span>
            </div>
            <div class="bg-white p-5 rounded-xl shadow-sm border border-slate-200/80">
                <span class="text-xs text-slate-400 font-semibold uppercase tracking-wider block">Fase II/III & Antiox</span>
                <span class="text-3xl font-bold text-amber-600 block mt-1" id="other-detox-count">0</span>
                <span class="text-xs text-slate-500 block mt-1">GSTs, UGTs, ABCs e SODs</span>
            </div>
        </div>

        <!-- CONTROLES -->
        <div class="bg-white p-6 rounded-xl shadow-sm border border-slate-200/80 mb-6">
            <div class="grid grid-cols-1 md:grid-cols-3 gap-4">
                <!-- Busca de Texto -->
                <div>
                    <label class="block text-xs font-semibold text-slate-500 uppercase tracking-wider mb-2">Procurar Termo ou ID</label>
                    <div class="relative">
                        <input type="text" id="search-input" placeholder="Ex: CYP6, GST, rna-gnl..." 
                               class="w-full bg-slate-50 border border-slate-300 rounded-lg px-3 py-2 pl-10 text-sm focus:outline-none focus:ring-2 focus:ring-emerald-500 focus:border-emerald-500">
                        <i class="fa-solid fa-magnifying-glass absolute left-3.5 top-3 text-slate-400 text-xs"></i>
                    </div>
                </div>
                <!-- Filtro de Família -->
                <div>
                    <label class="block text-xs font-semibold text-slate-500 uppercase tracking-wider mb-2">Filtrar por Família de Detox</label>
                    <select id="family-filter" class="w-full bg-slate-50 border border-slate-300 rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-emerald-500 focus:border-emerald-500">
                        <option value="all">Todas as Famílias</option>
                        <option value="Fase I: Cytochrome P450 (CYP)">Fase I: Cytochrome P450 (CYP)</option>
                        <option value="Fase I: Carboxylesterase (CCE)">Fase I: Carboxylesterase (CCE)</option>
                        <option value="Fase II: Glutathione S-transferase (GST)">Fase II: Glutathione S-transferase (GST)</option>
                        <option value="Fase II: UDP-glucuronosyltransferase (UGT)">Fase II: UDP-glucuronosyltransferase (UGT)</option>
                        <option value="Fase III: Transportadores ABC">Fase III: Transportadores ABC</option>
                        <option value="Antioxidante: Superoxide Sodumase (SOD)">Antioxidantes</option>
                    </select>
                </div>
                <!-- Filtro de On-Score -->
                <div>
                    <label class="block text-xs font-semibold text-slate-500 uppercase tracking-wider mb-2">Consistência (On-Score)</label>
                    <select id="score-filter" class="w-full bg-slate-50 border border-slate-300 rounded-lg px-3 py-2 text-sm focus:outline-none focus:ring-2 focus:ring-emerald-500 focus:border-emerald-500">
                        <option value="all">Qualquer Score (&ge; 2)</option>
                        <option value="3">Exatamente 3 (Resposta Core Universal)</option>
                        <option value="2">Exatamente 2 (Ativo em dois estudos)</option>
                    </select>
                </div>
            </div>
            
            <!-- Ações Rápidas de Cópia -->
            <div class="mt-4 pt-4 border-t border-slate-100 flex flex-wrap gap-2 justify-end">
                <button onclick="copiarLatexTabela()" class="bg-slate-800 text-slate-100 hover:bg-slate-900 px-4 py-2 rounded-lg text-xs font-semibold inline-flex items-center gap-1.5 transition-colors cursor-pointer">
                    <i class="fa-solid fa-code"></i> Copiar Código LaTeX (booktabs)
                </button>
                <button onclick="copiarMarkdownTabela()" class="bg-emerald-600 text-white hover:bg-emerald-700 px-4 py-2 rounded-lg text-xs font-semibold inline-flex items-center gap-1.5 transition-colors cursor-pointer">
                    <i class="fa-solid fa-paste"></i> Copiar Tabela Markdown
                </button>
            </div>
        </div>

        <!-- TABELA INTERATIVA -->
        <div class="bg-white rounded-xl shadow-sm border border-slate-200/80 overflow-hidden">
            <div class="overflow-x-auto">
                <table class="w-full text-left border-collapse">
                    <thead>
                        <tr class="bg-slate-100 border-b border-slate-200 text-slate-500 text-xs font-semibold uppercase tracking-wider">
                            <th class="py-3.5 px-4">Transcript ID</th>
                            <th class="py-3.5 px-4 w-[35%]">Anotação Funcional</th>
                            <th class="py-3.5 px-4">Família de Detox</th>
                            <th class="py-3.5 px-4 text-center">Abamectina<br><span class="text-[10px] font-medium lowercase">log2 ratio</span></th>
                            <th class="py-3.5 px-4 text-center">Spinosad<br><span class="text-[10px] font-medium lowercase">log2 ratio</span></th>
                            <th class="py-3.5 px-4 text-center">Emamectina<br><span class="text-[10px] font-medium lowercase">log2 ratio</span></th>
                            <th class="py-3.5 px-4 text-center">On-Score<br><span class="text-[10px] font-medium lowercase">UP</span></th>
                        </tr>
                    </thead>
                    <tbody id="table-body" class="divide-y divide-slate-100 text-sm">
                        <!-- Gerado dinamicamente via JS -->
                    </tbody>
                </table>
            </div>
            <!-- Estado Vazio -->
            <div id="no-results" class="hidden p-8 text-center text-slate-400">
                <i class="fa-solid fa-folder-open text-3xl mb-2"></i>
                <p>Nenhum transcrito encontrado com os filtros selecionados.</p>
            </div>
        </div>
    </main>

    <!-- NOTIFICAÇÃO -->
    <div id="toast" class="fixed bottom-5 right-5 bg-slate-900 text-white px-4 py-3 rounded-xl shadow-lg border border-slate-700 text-sm hidden items-center gap-2 transition-all">
        <i class="fa-solid fa-circle-check text-emerald-400"></i>
        <span id="toast-text">Tabela copiada para a área de transferência!</span>
    </div>

    <!-- SCRIPT DE INTERAÇÃO -->
    <script>
        // Dados injetados dinamicamente via Python
        const dadosOriginais = {dados_json};

        function renderizarTabela(filtros = {{}}) {{
            const tbody = document.getElementById('table-body');
            const noResults = document.getElementById('no-results');
            tbody.innerHTML = '';
            
            const txt = filtros.text ? filtros.text.toLowerCase() : '';
            const family = filtros.family || 'all';
            const score = filtros.score || 'all';
            
            let filtrados = dadosOriginais.filter(item => {{
                // Filtro de Texto
                const matchText = item.id.toLowerCase().includes(txt) || item.annotation.toLowerCase().includes(txt);
                // Filtro de Família
                const matchFamily = (family === 'all') || (item.family === family);
                // Filtro de Score
                const matchScore = (score === 'all') || (item.score === parseInt(score));
                
                return matchText && matchFamily && matchScore;
            }});

            // Ordenação por padrão: On-Score Decrescente, depois pelo maior Fold-Change
            filtrados.sort((a, b) => {{
                if (b.score !== a.score) return b.score - a.score;
                return Math.max(b.fc_abam, b.fc_spin, b.fc_emam) - Math.max(a.fc_abam, a.fc_spin, a.fc_emam);
            }});

            if (filtrados.length === 0) {{
                noResults.classList.remove('hidden');
            }} else {{
                noResults.classList.add('hidden');
                
                filtrados.forEach(item => {{
                    const tr = document.createElement('tr');
                    
                    // Destaca linhas com Score 3 (resposta mestre)
                    if (item.score === 3) {{
                        tr.className = 'bg-indigo-50/50 hover:bg-indigo-50 font-medium border-l-4 border-l-indigo-500';
                    }} else {{
                        tr.className = 'hover:bg-slate-50/80';
                    }}
                    
                    const badgeDetox = obterBadgeDetox(item.family);
                    const classAbam = item.fc_abam >= 1.0 ? 'text-emerald-600 font-semibold' : 'text-slate-500';
                    const classSpin = item.fc_spin >= 1.0 ? 'text-emerald-600 font-semibold' : 'text-slate-500';
                    const classEmam = item.fc_emam >= 1.0 ? 'text-emerald-600 font-semibold' : 'text-slate-500';
                    
                    tr.innerHTML = `
                        <td class="py-3 px-4 font-mono text-xs text-slate-600">${{item.id}}</td>
                        <td class="py-3 px-4 text-slate-800 break-words">${{item.annotation}}</td>
                        <td class="py-3 px-4">${{badgeDetox}}</td>
                        <td class="py-3 px-4 text-center ${{classAbam}}">${{item.fc_abam.toFixed(2)}}</td>
                        <td class="py-3 px-4 text-center ${{classSpin}}">${{item.fc_spin.toFixed(2)}}</td>
                        <td class="py-3 px-4 text-center ${{classEmam}}">${{item.fc_emam.toFixed(2)}}</td>
                        <td class="py-3 px-4 text-center">
                            <span class="inline-flex items-center justify-center px-2 py-1 rounded-full text-xs font-bold ${{item.score === 3 ? 'bg-indigo-600 text-white' : 'bg-slate-200 text-slate-700'}}">
                                ${{item.score}}
                            </span>
                        </td>
                    `;
                    tbody.appendChild(tr);
                }});
            }}
            
            // Atualiza cartões de métricas
            document.getElementById('total-count').textContent = filtrados.length;
            document.getElementById('core-count').textContent = filtrados.filter(x => x.score === 3).length;
            document.getElementById('phase1-count').textContent = filtrados.filter(x => x.family.includes('CYP') || x.family.includes('CCE')).length;
            document.getElementById('other-detox-count').textContent = filtrados.filter(x => x.family.includes('GST') || x.family.includes('UGT') || x.family.includes('ABC') || x.family.includes('SOD') || x.family.includes('CAT') || x.family.includes('Peroxidase')).length;
        }}

        function obterBadgeDetox(fam) {{
            if (fam.includes('CYP')) return '<span class="bg-purple-100 text-purple-800 text-xs px-2.5 py-0.5 rounded-full font-medium inline-block border border-purple-200"><i class="fa-solid fa-dna mr-1"></i> P450 (Fase I)</span>';
            if (fam.includes('CCE')) return '<span class="bg-blue-100 text-blue-800 text-xs px-2.5 py-0.5 rounded-full font-medium inline-block border border-blue-200"><i class="fa-solid fa-droplet mr-1"></i> CCE (Fase I)</span>';
            if (fam.includes('GST')) return '<span class="bg-emerald-100 text-emerald-800 text-xs px-2.5 py-0.5 rounded-full font-medium inline-block border border-emerald-200"><i class="fa-solid fa-shield mr-1"></i> GST (Fase II)</span>';
            if (fam.includes('UGT')) return '<span class="bg-amber-100 text-amber-800 text-xs px-2.5 py-0.5 rounded-full font-medium inline-block border border-amber-200"><i class="fa-solid fa-capsules mr-1"></i> UGT (Fase II)</span>';
            if (fam.includes('ABC')) return '<span class="bg-rose-100 text-rose-800 text-xs px-2.5 py-0.5 rounded-full font-medium inline-block border border-rose-200"><i class="fa-solid fa-truck-ramp-box mr-1"></i> ABC (Fase III)</span>';
            if (fam.includes('SOD') || fam.includes('CAT') || fam.includes('Peroxidase')) return '<span class="bg-teal-100 text-teal-800 text-xs px-2.5 py-0.5 rounded-full font-medium inline-block border border-teal-200"><i class="fa-solid fa-snowflake mr-1"></i> Antiox</span>';
            return '<span class="bg-slate-100 text-slate-500 text-xs px-2.5 py-0.5 rounded-full font-medium inline-block">Outros Processos</span>';
        }}

        function mostrarToast(texto) {{
            const toast = document.getElementById('toast');
            document.getElementById('toast-text').textContent = texto;
            toast.classList.remove('hidden');
            toast.classList.add('flex');
            setTimeout(() => {{
                toast.classList.remove('flex');
                toast.classList.add('hidden');
            }}, 3000);
        }}

        function copiarLatexTabela() {{
            let tex = `\\\\begin{{table}}[htbp]
  \\\\centering
  \\\\caption{{Tabela de metanálise filtrada por On-Score em Phthorimaea absoluta.}}
  \\\\label{{tab:metanalise_detox}}
  \\\\small
  \\\\begin{{tabular}}{{l p{{6cm}} c c c c c}}
    \\\\toprule
    \\\\textbf{{Transcript ID}} & \\\\textbf{{Anotação Funcional}} & \\\\textbf{{Abam.}} & \\\\textbf{{Spin.}} & \\\\textbf{{Emam.}} & \\\\textbf{{On-Score}} & \\\\textbf{{On-Ratio}} \\\\\\\\
    \\\\midrule\\n`;

            dadosOriginais.forEach(item => {{
                // Escapa underscores do ID e anotação para o LaTeX compilado
                const escId = item.id.replace(/_/g, '\\\\_');
                const escAnnot = item.annotation.replace(/_/g, '\\\\_').substring(0, 60) + (item.annotation.length > 60 ? '...' : '');
                
                tex += `    ${{escId}} & ${{escAnnot}} & ${{item.fc_abam.toFixed(2)}} & ${{item.fc_spin.toFixed(2)}} & ${{item.fc_emam.toFixed(2)}} & ${{item.score}} & ${{item.ratio.toFixed(2)}} \\\\\\\\\\n`;
            }});

            tex += `    \\\\bottomrule
  \\\\end{{tabular}}
\\\\end{{table}}`;

            const area = document.createElement('textarea');
            area.value = tex;
            document.body.appendChild(area);
            area.select();
            document.execCommand('copy');
            document.body.removeChild(area);
            mostrarToast("Código LaTeX profissional (booktabs) copiado!");
        }}

        function copiarMarkdownTabela() {{
            let md = "| Transcript ID | Anotação Funcional | Log2FC Abam. | Log2FC Spin. | Log2FC Emam. | On-Score (UP) | On-Ratio |\\n|---|---|---|---|---|---|---|\\n";
            dadosOriginais.forEach(item => {{
                md += `| \`${{item.id}}\` | ${{item.annotation}} | ${{item.fc_abam.toFixed(2)}} | ${{item.fc_spin.toFixed(2)}} | ${{item.fc_emam.toFixed(2)}} | **${{item.score}}** | ${{item.ratio.toFixed(2)}} |\\n`;
            }});
            
            const area = document.createElement('textarea');
            area.value = md;
            document.body.appendChild(area);
            area.select();
            document.execCommand('copy');
            document.body.removeChild(area);
            mostrarToast("Tabela em formato Markdown copiada!");
        }}

        // Eventos de Escuta
        document.getElementById('search-input').addEventListener('input', (e) => {{
            renderizarTabela({{
                text: e.target.value,
                family: document.getElementById('family-filter').value,
                score: document.getElementById('score-filter').value
            }});
        }});

        document.getElementById('family-filter').addEventListener('change', (e) => {{
            renderizarTabela({{
                text: document.getElementById('search-input').value,
                family: e.target.value,
                score: document.getElementById('score-filter').value
            }});
        }});

        document.getElementById('score-filter').addEventListener('change', (e) => {{
            renderizarTabela({{
                text: document.getElementById('search-input').value,
                family: document.getElementById('family-filter').value,
                score: e.target.value
            }});
        }});

        // Carregamento inicial da página
        renderizarTabela();
    </script>
</body>
</html>
"""
    
    with open(caminho_saida, 'w', encoding='utf-8') as f:
        f.write(html_content)
        
    print(f"   ✅ Dashboard interativo guardado em: '{caminho_saida}'")

def executar_gerador_tabelas():
    print("==================================================")
    print("      SAQE - COMPILADOR DE TABELAS ACADÉMICAS     ")
    print("==================================================")
    
    input_csv = "quant/metanalise_completa_anotada.csv"
    
    if not os.path.exists(input_csv):
        print(f"❌ [ERRO] O ficheiro anotado da metanálise '{input_csv}' não foi encontrado!")
        print("Dica: Corra primeiro o script de anotação funcional: 'python3 anotar_resultados.py'")
        return
        
    df = pd.read_csv(input_csv, index_col=0)
    total_linhas = len(df)
    
    # Filtra apenas os genes biologicamente relevantes (On-Score >= 2) para a tabela de publicação principal
    # Isso evita entupir o manuscrito com milhares de genes que flutuaram em apenas 1 comparação
    df_filtrado = df[(df['Bono_UP_Score'] >= 2) | (df['Bono_DOWN_Score'] >= 2)].copy()
    
    print(f"📊 Registos totais anotados: {total_linhas}")
    print(f"📊 Selecionados {len(df_filtrado)} candidatos consistentes (On-Score >= 2) para as tabelas principais.")
    
    # Cria pasta quant se necessário
    os.makedirs("quant", exist_ok=True)
    
    # Exportação dos formatos
    gerar_tabela_latex(df_filtrado, "quant/tabela_artigo_latex.tex")
    gerar_tabela_markdown(df_filtrado, "quant/tabela_artigo_markdown.md")
    gerar_dashboard_html(df_filtrado, "quant/tabela_interativa_artigo.html")
    
    print("\n==================================================")
    print("🎉 COMPILAÇÃO FINALIZADA COM SUCESSO!")
    print("==================================================")
    print("📂 Ficheiros gerados na pasta 'quant/':")
    print("   👉 'tabela_artigo_latex.tex'    -> Código LaTeX (Overleaf/Texmaker)")
    print("   👉 'tabela_artigo_markdown.md' -> Tabela para relatórios ou cópia MS Word")
    print("   👉 'tabela_interativa_artigo.html' -> Dashboard para explorar, pesquisar e filtrar")
    print("\nDica de ouro: Abra o ficheiro 'tabela_interativa_artigo.html' em qualquer navegador")
    print("para copiar códigos de partes da tabela filtrados em poucos segundos!")

if __name__ == "__main__":
    executar_gerador_tabelas()
