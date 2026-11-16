# Verificação curta da geometria e das raízes

Antes de integrar os 18 resultados no artigo, verificar duas possíveis falhas que mudariam a conclusão: erro no fator de superfície da parametrização e localização incorreta das fronteiras. Sem novos casos de desempenho ou ajuste de tolerância.

1. Com `q=0` e `Rx=Ry`, a secção reduz-se ao cilindro do teste 45. Usar as duas alturas já existentes que apresentaram falso seguro no cilindro; comparar a referência mpmath do novo código com o valor já guardado na referência circular. Critério: erro absoluto inferior a 1e−12 pN, muito abaixo da tolerância da força. Este controlo avalia fator, sinal e lei, não é nova validação independente da física.
2. Nas 18 alturas, comparar as raízes Brent com a resolução da equação `2q c²+Ry c+h−q−B=0`, c=cosθ, para as raízes fisicamente admissíveis em (−1,1), B∈{0,a}. Reportar o maior desvio angular e eventuais discrepâncias do número de raízes. A forma algébrica existe para esta geometria; a busca numérica foi escolhida para medir seu custo e não por impossibilidade analítica.

Parar após estes controlos leves, CPU1. A comparação de custos continua a ser a grade previamente fixada, sem novos métodos de referência nem redefinição de tolerâncias.
