from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from scipy.stats import t

st.set_page_config(page_title='PE / Graphene · MD Explorer', page_icon='⬡', layout='wide')
ROOT=Path(__file__).parent
DATA=ROOT/'data'
TEAL='#1FB5C4'; AMBER='#F6AD55'; COLORS=[TEAL,AMBER,'#9C9BFF']
st.markdown('''<style>
html,body,[data-testid="stApp"]{font-family:'Segoe UI','Open Sans',sans-serif}
.block-container{max-width:1480px;padding:4rem 3rem 3rem}
[data-testid="stAppViewContainer"]{background:#041522}
[data-testid="stMarkdownContainer"] p,[data-testid="stMarkdownContainer"] li{font-size:19px;line-height:1.7}
h1{font-size:clamp(36px,3.5vw,56px)!important;letter-spacing:-1.8px;font-weight:700!important;line-height:1.12!important;margin-bottom:22px!important}
h2{font-size:30px!important;letter-spacing:-.7px}h3{font-size:25px!important}
[data-testid="stCaptionContainer"] p{font-size:15px!important;line-height:1.6;color:#a3bac9}
[data-testid="stSidebar"]{background:#092335;border-right:1px solid #224051}
[data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p{font-size:17px}
[data-testid="stSidebar"] [role="radiogroup"] label{padding:7px 10px;border-radius:7px;transition:background .18s ease,transform .18s ease}
[data-testid="stSidebar"] [role="radiogroup"] label:hover{background:#12364b;transform:translateX(3px)}
[data-testid="stMetric"]{background:#0b2638;border:1px solid #234555;border-top:2px solid #1FB5C4;border-radius:10px;padding:23px;transition:transform .2s ease,border-color .2s ease;animation:arrive .55s ease both}
[data-testid="stMetric"]:hover{transform:translateY(-4px);border-color:#1FB5C4}
[data-testid="stMetricLabel"] p{font-size:17px!important}
[data-testid="stMetricValue"]{font-size:46px!important;font-weight:600}
.hero{position:relative;overflow:hidden;padding:30px 34px;border-left:3px solid #1FB5C4;background:linear-gradient(110deg,#102d40,#071c2c);border-radius:0 12px 12px 0;margin-bottom:26px;animation:arrive .65s ease both}
.hero h2{font-size:clamp(24px,2.4vw,34px)!important;line-height:1.3;margin:0 0 16px}
.hero p{max-width:850px;margin-bottom:14px}.hero span{font-size:16px}
.hero::after{content:"";position:absolute;right:-55px;top:-80px;width:220px;height:220px;border:1px solid #1fb5c428;border-radius:50%;pointer-events:none}
.note{border-left:3px solid #F6AD55;padding:18px 22px;background:#0b2638;border-radius:0 8px 8px 0;margin:22px 0;font-size:18px;line-height:1.7}
.stButton button,.stDownloadButton button{font-size:16px;border-radius:8px;transition:transform .2s ease,border-color .2s ease}
.stButton button:hover,.stDownloadButton button:hover{transform:translateY(-2px);border-color:#1FB5C4}
[data-testid="stImage"] img{border-radius:10px}
@keyframes arrive{from{opacity:0;transform:translateY(12px)}to{opacity:1;transform:translateY(0)}}
@media(prefers-reduced-motion:reduce){*,*::before,*::after{animation:none!important;transition:none!important;transform:none!important}}
@media(max-width:760px){.block-container{padding:3.5rem 1rem 2rem}.hero{padding:23px}h1{font-size:36px!important}}
</style>''',unsafe_allow_html=True)

def read(name):
    try: return pd.read_csv(DATA/name)
    except Exception as e: st.error(f'Cannot read {name}: {e}'); st.stop()

def results(name):
    d=read(name)
    required=['run','n_flakes','wt_percent','Tg_K','E_avg','Y_avg']
    if not set(required)<=set(d): st.error(f'{name}: required columns {required}'); st.stop()
    if d.empty or d['run'].duplicated().any(): st.error(f'{name}: empty data or duplicate run IDs.'); st.stop()
    for col in required[1:]:
        d[col]=pd.to_numeric(d[col],errors='coerce')
    if not np.isfinite(d[required[1:]].to_numpy()).all(): st.error(f'{name}: missing or nonnumeric measurements.'); st.stop()
    return d

def draw(fig):
    fig.update_layout(template='plotly_dark',paper_bgcolor='#031323',plot_bgcolor='#031323',font=dict(color='#E8F4FA',family='Segoe UI',size=17),colorway=COLORS,margin=dict(l=25,r=25,t=45,b=40),height=440,legend=dict(orientation='h',y=1.12),hovermode='closest')
    fig.update_xaxes(gridcolor='#183248',zeroline=False)
    fig.update_yaxes(gridcolor='#183248',zeroline=False)
    st.plotly_chart(fig,use_container_width=True,theme=None)

def stats(d,col):
    s=d.groupby('wt_percent')[col].agg(['mean','std','count']).reset_index()
    s['CI95']=s.apply(lambda r: t.ppf(.975,r['count']-1)*r['std']/np.sqrt(r['count']) if r['count']>1 else np.nan,axis=1)
    return s

def property_chart(d,col,label):
    s=stats(d,col);f=go.Figure()
    f.add_trace(go.Scatter(x=d.wt_percent,y=d[col],mode='markers',name='Individual seed',text=d.run,marker=dict(color=AMBER,size=7,opacity=.65),hovertemplate='%{text}<br>%{x:.3f} wt%<br>%{y:.4f}<extra></extra>'))
    f.add_trace(go.Scatter(x=s.wt_percent,y=s['mean'],error_y=dict(type='data',array=s.CI95),mode='lines+markers',name='Mean · 95% CI',line=dict(color=TEAL,width=3),marker=dict(size=11)))
    f.update_layout(xaxis_title='Graphene content (wt%)',yaxis_title=label)
    draw(f);return s

small=results('results_xyz.csv');big=results('results_big.csv')
st.sidebar.markdown('## ⬡ PE / GRAPHENE')
st.sidebar.caption('Molecular dynamics study')
page=st.sidebar.radio('Explore',['Overview','Systems & Method','Stress–Strain','Modulus / Yield / Tg','Optimum filler content','Interface','Theory · Halpin–Tsai','Team'],key='page')
size=st.sidebar.selectbox('Study',['Small system','Large system'])
d=small if size=='Small system' else big
st.sidebar.caption(f'{len(d)} configurations · {d.n_flakes.nunique()} compositions')
if st.sidebar.button('Reload data'): st.rerun()
st.sidebar.download_button('Download active CSV',d.to_csv(index=False),'active_results.csv','text/csv')
st.sidebar.caption('CSV values are read on every rerun. Directions are averaged within each seed; they are not independent replicas.')
st.title(page if page!='Overview' else 'Polyethylene + graphene')

if page=='Overview':
    st.markdown("""<style>
    .home-intro{font-size:22px;color:#b5cbd8;max-width:820px;line-height:1.6;margin-bottom:24px}
    .live-tag{display:inline-flex;align-items:center;gap:9px;background:#103c47;border:1px solid #226373;border-radius:30px;padding:7px 14px;font-size:13px;letter-spacing:1px;color:#84e6ec;margin-bottom:18px}
    .live-dot{height:7px;width:7px;border-radius:50%;background:#4dd7df;animation:glow 2.5s ease-in-out infinite}
    @keyframes glow{50%{box-shadow:0 0 0 6px #1fb5c422}}
    .result-card{padding:20px 22px;background:linear-gradient(125deg,#123649,#0b2334);border:1px solid #285168;border-radius:12px;min-height:160px;animation:arrive .45s ease both}
    .result-label{font-size:15px;color:#b1cbd9;margin-bottom:9px}
    .result-number{font-size:36px;font-weight:650;letter-spacing:-1px;color:#f1f9fc}
    .result-unit{font-size:17px;color:#a5bfce}
    .result-change{color:#f6ad55;font-size:16px;margin-top:8px}
    </style>""",unsafe_allow_html=True)
    st.markdown('<div class="live-tag"><span class="live-dot"></span>EXPLORE THE SIMULATION RESULTS</div><div class="home-intro">What changes when graphene enters the polymer?<br>Choose a composition and compare its response with neat polyethylene.</div>',unsafe_allow_html=True)
    left,right=st.columns([1.15,1],gap='large')
    with left:
        st.subheader('Choose your material')
        options=sorted(d.n_flakes.unique())
        flakes=st.select_slider('Graphene flakes',options=options,value=options[-1],format_func=lambda n:'Neat PE' if n==0 else f'{int(n)} flake'+('s' if n>1 else ''))
        subset=d[d.n_flakes==flakes]
        base=d[d.n_flakes==0]
        wt=subset.wt_percent.mean()
        st.markdown(f'**{size} · {wt:.3f} wt% graphene**')
        st.caption(f'{len(subset)} independent configurations in this composition. Change the study in the sidebar to compare model scales.')
        st.markdown('**Inspect the response**')
        feature=st.radio('Property',['Stiffness','Yield stress','Glass transition'],horizontal=True,label_visibility='collapsed')
        field={'Stiffness':'E_avg','Yield stress':'Y_avg','Glass transition':'Tg_K'}[feature]
        st.caption('Every value comes from the uploaded CSV. Directional mechanical results are averaged within each configuration.')
    with right:
        st.image(str(DATA/'big_2flakes.png'),use_container_width=True)
        st.caption('Reference structure: large system, two flakes. This image stays fixed as you change composition; it is not a live molecular simulation.')
    cards=st.columns(3)
    for box,col,title,unit in zip(cards,['E_avg','Y_avg','Tg_K'],['Young’s modulus','Yield stress','Glass transition'],['GPa','GPa','K']):
        val=subset[col].mean();ref=base[col].mean()
        delta=val-ref
        change=(f'{delta:+.1f} K vs neat PE' if col=='Tg_K' else f'{100*delta/ref:+.1f}% vs neat PE') if len(base) and ref!=0 else 'Baseline unavailable'
        if flakes==0: change='Neat-polymer reference'
        digits=1 if col=='Tg_K' else 3
        with box: st.markdown(f'<div class="result-card"><div class="result-label">{title}</div><div class="result-number">{val:.{digits}f} <span class="result-unit">{unit}</span></div><div class="result-change">{change}</div></div>',unsafe_allow_html=True)
    st.caption('Composition means. Changes shown above are descriptive, not claims of statistical significance.')
    st.subheader('How does it compare?')
    # Keep the independent seed observations and uncertainty visible alongside the means.
    compare=d[d.n_flakes.isin([0,flakes])]
    property_chart(compare,field,{'E_avg':'Young’s modulus (GPa)','Y_avg':'Yield stress (GPa)','Tg_K':'Glass transition (K)'}[field])
    st.caption('Dots: independent configurations. Whiskers: 95% Student-t confidence intervals of each mean. Selecting neat PE shows the reference alone.')
    def go_page(target): st.session_state.page=target
    st.subheader('Go deeper')
    x,y,z=st.columns(3)
    x.button('Explore tensile curves →',on_click=go_page,args=('Stress–Strain',),use_container_width=True)
    y.button('Compare all compositions →',on_click=go_page,args=('Modulus / Yield / Tg',),use_container_width=True)
    z.button('Inspect the interface →',on_click=go_page,args=('Interface',),use_container_width=True)
    st.markdown('<div class="note">A favourable polymer–graphene interface does not automatically produce measurable reinforcement. Explore both the response and its uncertainty.</div>',unsafe_allow_html=True)
elif page=='Systems & Method':
    st.subheader('One matrix. Two model scales.')
    st.dataframe(pd.DataFrame({'Study':['Small','Large'],'Polymer beads':[1000,4000],'Beads per chain':[50,50],'Atoms per flake':[42,150],'Configurations in CSV':[len(small),len(big)]}),hide_index=True,use_container_width=True)
    st.write('Build polymer + filler → minimize → NVT/NPT equilibration → cooling → independent x/y/z tensile tests → per-seed property extraction.')
    st.info('Protocol from the supplied project description: amorphous united-atom polyethylene; tensile testing at 100 K, nominal strain rate 10⁹ s⁻¹. These CSVs alone do not validate force fields or equilibration.')
    st.warning('The larger study changes both box size and flake size. It is not an isolated box-size convergence test.')
    st.dataframe(d.groupby(['n_flakes','wt_percent']).size().rename('Completed seeds').reset_index(),hide_index=True)
elif page=='Stress–Strain':
    st.caption('Provided representative large-system y-direction curves. This page is independent of the study selector; these are not ensemble averages.')
    selected=st.multiselect('Curves',['Neat PE','2 flakes'],default=['Neat PE','2 flakes'])
    window=st.slider('Display-only moving average (points)',1,201,31,step=2)
    raw=st.checkbox('Show raw measurements',True);f=go.Figure()
    for i,label in enumerate(selected):
        name='neat_ss_y.txt' if label=='Neat PE' else 'flakes2_ss_y.txt'
        q=pd.read_csv(DATA/name,sep=r'\s+',comment='#',header=None,names=['step','strain','stress','density'])
        if raw: f.add_trace(go.Scatter(x=q.strain,y=q.stress,name=label+' raw',line=dict(color=COLORS[i],width=1),opacity=.2))
        f.add_trace(go.Scatter(x=q.strain,y=q.stress.rolling(window,center=True,min_periods=1).mean(),name=label,line=dict(color=COLORS[i],width=2.5)))
    f.update_layout(xaxis_title='Engineering strain',yaxis_title='Reported tensile stress (GPa)');draw(f)
    st.caption('Smoothing changes the display only. Modulus and yield values elsewhere come from the supplied result CSVs, not this smoothing control.')
elif page=='Modulus / Yield / Tg':
    label=st.selectbox('Property',['Young’s modulus (GPa)','Yield stress (GPa)','Glass transition (K)'])
    col=dict(zip(['Young’s modulus (GPa)','Yield stress (GPa)','Glass transition (K)'],['E_avg','Y_avg','Tg_K']))[label]
    s=property_chart(d,col,label)
    st.caption('Two-sided Student-t confidence intervals for each composition mean: mean ± t₀.₉₇₅,n−1 × SD/√n. These are not confidence intervals of differences. No CI is available for n=1.')
    st.dataframe(s,hide_index=True,use_container_width=True)
    with st.expander('Inspect individual configurations'): st.dataframe(d,hide_index=True,use_container_width=True)
elif page=='Optimum filler content':
    st.subheader('A highest mean is not an established optimum')
    col=st.selectbox('Ranking property',['E_avg','Y_avg','Tg_K'])
    s=property_chart(d,col,{'E_avg':'Modulus (GPa)','Y_avg':'Yield stress (GPa)','Tg_K':'Tg (K)'}[col])
    best=s.loc[s['mean'].idxmax()]
    st.metric('Highest observed mean at',f'{best.wt_percent:.3f} wt%')
    st.info('Descriptive ranking only. An optimum requires a defined design objective, uncertainty analysis and adequate sampling of filler content. No validated optimum is established by this dashboard.')
    base=s.loc[s.wt_percent==0,'mean']
    if len(base):
        s['Difference vs neat']=s['mean']-base.iloc[0]
        if col!='Tg_K': s['Change (%)']=100*s['Difference vs neat']/base.iloc[0]
    st.dataframe(s,hide_index=True)
elif page=='Interface':
    st.caption('Small-system interface data supplied separately; not updated by the large-system selector.')
    q=read('interface_profile.csv');f=go.Figure()
    for i in range(1,4):
        mean=q[f'g{i}_mean'];sd=q[f'g{i}_sd'];color=COLORS[i-1]
        f.add_trace(go.Scatter(x=q.r_A,y=mean+sd,line=dict(width=0),showlegend=False,hoverinfo='skip'))
        f.add_trace(go.Scatter(x=q.r_A,y=mean-sd,fill='tonexty',fillcolor=['rgba(31,181,196,.15)','rgba(246,173,85,.15)','rgba(156,155,255,.15)'][i-1],line=dict(width=0),showlegend=False,hoverinfo='skip'))
        f.add_trace(go.Scatter(x=q.r_A,y=mean,name=f'{i} flake(s)',line=dict(color=color)))
    f.add_hline(y=1,line_dash='dot');f.update_layout(xaxis_title='Distance r (Å)',yaxis_title='Normalized polymer density');draw(f)
    st.caption('Bands are the supplied ±1 SD, not 95% CIs. Distance definition and normalization follow the source analysis.')
    a=read('adhesion_all.csv');f=go.Figure()
    f.add_trace(go.Box(x=a.n_flakes,y=a.E_per_flake_kcal_mol,boxpoints='all',jitter=.2,marker_color=TEAL,name='Per-flake interaction'))
    f.update_layout(xaxis_title='Number of flakes',yaxis_title='Interaction energy per flake (kcal/mol)');draw(f)
    st.info('Negative interaction energy indicates favourable interaction under the chosen model. It is not a direct measurement of interfacial strength.')
    st.dataframe(a,hide_index=True)
elif page=='Theory · Halpin–Tsai':
    st.subheader('Continuum model, explicit assumptions')
    st.latex(r'E_c/E_m=(1+\xi\eta V_f)/(1-\eta V_f),\quad \eta=(E_f/E_m-1)/(E_f/E_m+\xi)')
    st.write('A model comparison requires filler modulus, geometry/orientation factor and volume fraction. Weight fraction alone is insufficient.')
    st.info('Enter documented parameters below. No theory curve is presented as a verified prediction without those inputs.')
    a,b,c=st.columns(3)
    ef=a.number_input('Filler modulus (GPa); 0 = unset',min_value=0.,value=0.)
    xi=b.number_input('Geometry factor ξ; 0 = unset',min_value=0.,value=0.)
    rp=c.number_input('Polymer density (g/cm³); 0 = unset',min_value=0.,value=0.)
    rf=st.number_input('Filler density (g/cm³); 0 = unset',min_value=0.,value=0.)
    if min(ef,xi,rp,rf)>0:
        em=d.loc[d.wt_percent==0,'E_avg'].mean()
        w=np.linspace(0,d.wt_percent.max()/100,100);v=(w/rf)/(w/rf+(1-w)/rp);eta=(ef/em-1)/(ef/em+xi)
        f=go.Figure(go.Scatter(x=w*100,y=em*(1+xi*eta*v)/(1-eta*v),name='User-parameter theory'))
        s=stats(d,'E_avg');f.add_trace(go.Scatter(x=s.wt_percent,y=s['mean'],error_y=dict(array=s.CI95),mode='markers',name='MD mean · 95% CI'))
        f.update_layout(xaxis_title='Graphene (wt%)',yaxis_title='Modulus (GPa)');draw(f)
        st.caption('Exploratory model. Finite flakes, imperfect transfer and orientation can limit applicability. Parameter inputs are not fitted or validated here.')
elif page=='Team':
    st.subheader('Group 1 · IIT Indore')
    st.dataframe(pd.DataFrame({'Member':['Aashima','Abhinav','Akshat','Bharti'],'Proposal responsibility':['Structure building & force-field selection','Equilibration & deformation simulations','Data analysis & modulus extraction','Visualization & report/presentation']}),hide_index=True,use_container_width=True)
    st.caption('Responsibilities reproduce the original proposal; update them to reflect actual contributions before submission.')
st.divider();st.caption('PE / GRAPHENE · Molecular dynamics study')
