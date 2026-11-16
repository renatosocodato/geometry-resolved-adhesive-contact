# Teste previamente fixado: indicador geométrico versus Gauss–Kronrod adaptativo

16/11/2026. Objetivo: decidir se uma regra geométrica simples, no contacto de uma parede cilíndrica com um plano, acrescenta valor prático específico ao nicho de mecânica. O resultado anterior que origina a hipótese foi um falso seguro de comparação entre ordens globais perto de uma transição angular; não se assume que aconteça com Gauss–Kronrod adaptativo. Sem ajustar os casos ou indicador depois da execução.

## Problema e indicadores

Para raio R e comprimento L=0,75 µm, gap `g(θ)=h+R cos θ`, θ∈[0,π], a força normal da parede é `F=−2RL∫ψ'(g)dθ`. k=20 pN/µm³, w=0,002 pN/µm, a=0,1 µm e a lei por ramos é a lei existente do projeto; não há tampas. A transição B=0 ou B=a cruza a secção quando `c=(B−h)/R ∈ (−1,1)`, com raiz `θ_B=acos(c)`. Proximidade da borda angular `d_B=min(θ_B,π−θ_B)`. Indicador fixo: há risco se alguma transição cruzar e `d_B < π/16`. O valor corresponde a uma escala angular de resolução escolhida antes dos resultados; não é cota de erro nem limiar ótimo. A análise usa dimensões modeladas, não dados biológicos.

1. **GK:** `scipy.integrate.quad` adaptativo sem pontos de quebra.
2. **Indicador + GK:** localizar as raízes por `acos`, avaliar a regra acima e, se houver risco, fornecer todas as raízes interiores a `quad`; caso contrário aplicar GK sem pontos.
3. **Todas as raízes + GK:** fornecer todas as raízes interiores sempre. É o controlo de segmentação conhecido, não um método reivindicado como novo.

As três condições usam a mesma tolerância de força absoluta 10⁻⁸ pN, `epsrel=0`, `limit=100`; converter epsabs para a integral angular pelo fator `2RL`. Contar avaliações da lei, tempo mediano de três execuções, tempo global medido dentro de cada método incluindo `acos`, decisões e chamada de quad; uma execução de aquecimento sem registo. Mesmo processo, CPU1/BLAS1. Relatar contagem de falha real (erro absoluto >τ), falso seguro (estimativa GK ≤τ mas erro real >τ), riscos marcados, avaliação e custo. Medir diretamente todo custo de localização e controlo no caso analítico; não extrapolar custo para raízes em superfícies gerais. Tempos são exploratórios e sujeitos a carga da máquina.

## Grade nova e referência

R∈{0,28209479177387814; 0,36} µm. Para cada R, B∈{0;0,1}, extremo angular θ≈0 ou π, e distância angular ε∈{0,01;0,03;0,08;0,15;0,23;0,38} rad, escolher respetivamente `h=B−R cos ε` e `h=B+R cos ε`: 48 casos, incluindo compressão e adesão. Acrescentar para cada raio duas alturas `h=0,05` e `h=R+0,12`, total 52. Os pontos são geometria sintética separada dos usados no diagnóstico 31; alguns regimes foram escolhidos porque os resultados anteriores sugeriram uma possível vulnerabilidade, pelo que a grade é orientada e não uma amostra representativa.

Referência independente: integral por ramos via primitivas trigonométricas em mpmath com 60 dígitos. Verificar duas alturas pré-fixadas, uma de cada raio, por integração numérica mpmath direta da mesma lei em 45 dígitos: `B=0, extremo π, ε=0,03` e `B=0,1, extremo 0, ε=0,15`. A referência não é uma cota intervalar certificada. Um controlo adicional compara cada saída e o erro interno retornado por GK; não usar apenas a estimativa da rotina como verdade.

## Decisão

Um contributo de nicho é sustentado se o indicador demonstrar uma vantagem concreta (evitar falhas de GK ou reduzir custo para precisão semelhante) na família fixada e se houver uma mensagem não redundante perante segmentação conhecida. Um resultado sem vantagem não exclui contribuição conceptual: avaliar se a estratégia que melhor resolveu o erro ou a fronteira de aplicabilidade, com comparação honesta e uma implicação mecânica clara, constitui nota útil no nicho. Nunca afirmar superioridade geral, novidade absoluta, garantia de erro, custo em corpos 3D ou validação nativa. Se a evidência bastar para uma mensagem delimitada, redigir um artigo curto com antecedente mais próximo e incertezas; caso contrário, preservar um esboço de mensagem e próximos passos. Não selecionar pós hoc subconjunto para alegação confirmatória. Uma única grade, sem reexecutar após ver resultados, salvo correção de erro de implementação documentada e reexecução completa. Sem submissão automática.
