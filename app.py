import streamlit as st
from backend import *
import plotly.figure_factory as ff

# -- Set page config
apptitle = 'Simulador Aerofólio'

st.set_page_config(page_title=apptitle, page_icon="🍃")

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
naca_points = st.sidebar.slider('Pontos na superfície do aerofólio', 100, 400, 200)  # min, max, default
resolution = st.sidebar.slider('Resolução da malha', 50, 400, 50)
length = st.sidebar.slider('Comprimento do Canal (m)', 1.0, 8.0, 4.0)
height = st.sidebar.slider('Altura do canal (m)', 0.5, 4.0, 1.0)
flow_speed = st.sidebar.slider('Velocidade do escoamento do ar (m/s)', 5, 50, 8)
speed_condition = float(flow_speed)#*float(height)

alpha = st.sidebar.slider('Ângulo de ataque (Graus)', 0, 360, 5)

# st.sidebar.markdown('#### Whitened and band-passed data')
# whiten = st.sidebar.checkbox('Whiten?', value=True)
# freqrange = st.sidebar.slider('Band-pass frequency range (Hz)', min_value=10, max_value=2000, value=(30,400))


# # -- Create sidebar for Q-transform controls
# st.sidebar.markdown('#### Q-tranform plot')
# vmax = st.sidebar.slider('Colorbar Max Energy', 10, 500, 25)  # min, max, default
# qcenter = st.sidebar.slider('Q-value', 5, 120, 5)  # min, max, default
# qrange = (int(qcenter*0.8), int(qcenter*1.2))

if naca_foil != "":
    try:
        naca_study = BladeMesh(profNaca=naca_foil, n_of_points_naca=naca_points, length=length, diameter=height,bc_upperwall_value=speed_condition, alpha=alpha)
        naca_study.create_mesh(mesh_resolution=resolution)
        st.write('Malha Criada com Sucesso')
        naca_study.create_formulation()
        st.write('Formulação criada com sucesso')
        naca_study.compute_solution()
        st.write('Solução computada com sucesso')
        fig1, fig2 = naca_study.visualize_solutions('web')
        # naca_study.visualize_solutions(plot_type='matplotlib')
        with st.beta_expander("Visualizar Soluções Graficamente"):
            st.write(fig1)
            st.write(fig2)
        # mesh_coords = naca_study.mesh.coordinates()
        # naca_study.vx.vector().vec().array
        delta = 0.6
        npoints = int(delta*2*100)
        x = np.linspace(naca_study.baricenter[0] - delta, naca_study.baricenter[0] + delta, npoints)
        y = np.linspace(naca_study.baricenter[1] - delta, naca_study.baricenter[1] + delta, npoints)
        x0 = interpolate(Expression("x[0]"),V)
        x1 = interpolate(Expression("x[1]"),V)
        fig = ff.create_quiver(x, y, [naca_study.vx(i,j) for i in x for j in y], [naca_study.vy(i,j) for i in x for j in y],
                       scale=1e-2,
                       arrow_scale=.4,
                       name='quiver',
                       line_width=1)
        st.plotly_chart(fig, use_container_width=True)

    except Exception as e:
        if naca_foil == '':
            st.header('Por favor, inserir dados do aerofólio')
        else:
            st.write('Deu ruim sô sepá nem existe esse naca ai', e)
else:
    st.header('Por favor, inserir dados do aerofólio')