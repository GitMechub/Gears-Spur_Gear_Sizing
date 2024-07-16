import streamlit as st
for k, v in st.session_state.items():
    st.session_state[k] = v

from PIL import Image
import os
path = os.path.dirname(__file__)
my_file = path+'/images/mechub_logo.png'
img = Image.open(my_file)

st.set_page_config(
    page_title='References - Spur Gear Sizing',
    layout="wide",
    page_icon=img
                   )

st.sidebar.image(img)
st.sidebar.markdown("[![YouTube](https://img.shields.io/badge/YouTube-FF0000?style=for-the-badge&logo=youtube&logoColor=white)](https://www.youtube.com/@Mechub?sub_confirmation=1) [![GitHub](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white)](https://github.com/GitMechub)")

hide_menu = '''
        <style>
        #MainMenu {visibility: hidden; }
        footer {visibility: hidden;}
        </style>
        '''
st.markdown(hide_menu, unsafe_allow_html=True)

st.header("References", anchor=False, divider='gray')
st.write("1.   Norton, R. L. (2010). Machine Design: An Integrated Approach (4th ed). Prentice Hall.")
st.write("2.   MatWeb. (n.d.). MatWeb Material Property Data. Retrieved from http://www.matweb.com")
st.write("3.   ASM Handbook, Volume 19 - Fatigue and Fracture. ASM International Handbook Committee. ASM International. ISBN: 9781615031672. Available at: https://books.google.com.br/books/about/ASM_Handbook_Volume_19_Fatigue_and_Fract.html?id=FlwvygEACAAJ&redir_esc=y")
st.write("4.   Zou, Ting & Shaker, Mathew & Angeles, Jorge & Morozov, Alexei. (2014). Optimization of Tooth Root Profile of Spur Gears for Maximum Load-Carrying Capacity. 10.13140/2.1.4240.1286.")
st.write("5.   American Gear Manufacturers Association (AGMA)")
