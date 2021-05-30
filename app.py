import streamlit as st
from backend import *
import plotly.figure_factory as ff

# -- Set page config
apptitle = 'Simulador Aerofólio - MECÂNICA UFRJ'

st.set_page_config(page_title=apptitle, page_icon="🍃")

images = ['./logos/Marca_Principal.png','./logos/logo_dem.jpg']

cols = st.beta_columns(4)
cols[0].image(images[0], use_column_width=True)
cols[1].image(images[1], width=150)
# st.image(images, height=200)

# Title the app
st.title('Simulador Aerofólio')

st.markdown("""
 * Use o menu para preencher os parâmetros
 * A  análise será apresentada abaixo
""")

# -- Create Sidebar 
st.sidebar.markdown("## Selecione os dados para a simulação do aerofólio")
naca_foil = st.sidebar.text_input('Qual o perfil NACA será simulado?')

# -- Create sidebar for plot controls
st.sidebar.markdown('## Parâmetros do aerofólio')
naca_points = st.sidebar.slider('Pontos na superfície do aerofólio', 10, 400, 200)  # min, max, default
resolution = st.sidebar.slider('Resolução da malha', 50, 400, 50)
length = st.sidebar.slider('Comprimento da tubulação (m)', 1.0, 8.0, 2.0)
height = st.sidebar.slider('Altura do canal (m)', 0.5, 4.0, 1.2)
flow_speed = st.sidebar.slider('Velocidade do escoamento do ar (m/s)', 10, 400, 200)
circ_cc = st.sidebar.slider('Condição de Contorno no Aerofólio', 10, 400, int(flow_speed/2))
speed_condition = float(flow_speed)#*float(height)
alpha = st.sidebar.slider('Ângulo de ataque (Graus)', 0, 360, 5)
rho = st.sidebar.slider('Massa Específica do Ar (SI)', 0.8, 3.0, 1.292)

if naca_foil != "":
    try:
        naca_study = BladeMesh(profNaca=naca_foil, n_of_points_naca=naca_points, length=length, diameter=height,flow_speed=speed_condition, alpha=alpha, rho=rho)
        naca_study.create_mesh(mesh_resolution=resolution)
        naca_study.create_formulation(circ_cc=circ_cc)
        naca_study.compute_solution()
        fig1, fig2, fig3, fig4 = naca_study.visualize_solutions('web')
        circ, lift, l1 = naca_study.calculate_lift()
        # naca_study.visualize_solutions(plot_type='matplotlib')
        with st.beta_expander("Visualizar Soluções Graficamente"):
            st.write(fig1)
            st.write(fig2)
            st.write(fig3)
            st.write(fig4)
        # print(mpld3.fig_to_html(fig1).split('</style>\n\n')[-1][0:200])
        st.header('Teorema da força de sustentação:')
        st.markdown('A força de sustentação $L$ é dada por:')
        st.latex(r'''L = (\rho V)_{\infty} \Gamma ''')
        st.markdown('Em que $\Gamma$ é a circulação no aerofólio, dada a condição de Kutta:')
        st.latex(r"\Gamma =  \oint_{\text{Aerofólio}} \vec{V}\cdot d\vec{s}")
        st.markdown('Para os dados definidos, tem-se:')
        st.markdown(f'$\Gamma = {circ}$')
        st.markdown(f'$L = {lift}N/m$ ')
    except Exception as e:
        if naca_foil == '':
            st.header('Por favor, inserir dados do aerofólio')
        else:
            st.write('Por favor, verifique se esse aerofólio é NACA e pertence à serie 4 ou 5.', e)
else:
    st.header('Por favor, inserir dados do aerofólio')
    
