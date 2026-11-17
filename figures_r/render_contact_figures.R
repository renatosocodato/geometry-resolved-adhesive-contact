#!/usr/bin/env Rscript
# Vector-first representation of the existing analytical and benchmark results.
options(warn=1)
file_arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)[1]
here <- dirname(normalizePath(sub('^--file=', '', file_arg)))
paper <- dirname(here); input <- paper
suppressPackageStartupMessages({library(ggplot2);library(cowplot);library(svglite)})
cowplot::set_null_device('agg')
navy <- '#315778'; gold <- '#C29522'; grey <- '#888F97'; blue <- '#6D8DA6'; ink <- '#202B37'
pal <- c('Adaptive GK'=navy,'Guarded GK'=grey,'Split GK'=gold,'Closed form'=blue)
theme_clean <- function() theme_classic(base_size=10,base_family='Arial') + theme(
 text=element_text(colour=ink),axis.text=element_text(colour=ink,size=9.5),
 axis.title=element_text(size=10),axis.line=element_line(linewidth=.25),
 axis.ticks=element_line(linewidth=.25),legend.title=element_blank(),
 legend.text=element_text(size=9),legend.key.height=grid::unit(3,'mm'),
 legend.position='bottom',plot.margin=margin(5,5,5,5))
tau<-1e-8; R<-.28209479177387814; w<-.002; a<-.1; k<-20; h<-.32

# Figure 1a: vector cross-section; coordinate follows g(theta)=h+R cos(theta).
ang<-seq(0,2*pi,length.out=361)
circle<-data.frame(x=R*sin(ang),y=h+R*cos(ang))
theta<-2.27; px<-R*sin(theta); py<-h+R*cos(theta)
pa<-ggplot()+annotate('rect',xmin=-.48,xmax=.43,ymin=-.035,ymax=0,fill='#EFF1F4')+
 geom_hline(yintercept=0,linewidth=.5,colour=ink)+
 geom_path(data=circle,aes(x,y),linewidth=.65,colour=navy)+
 annotate('point',x=0,y=h,size=1.7,colour=ink)+
 annotate('segment',x=0,y=h,xend=px,yend=py,linewidth=.38,colour=ink,
          arrow=grid::arrow(length=grid::unit(1.5,'mm'),type='closed'))+
 annotate('point',x=px,y=py,size=2.1,colour=gold)+
 annotate('segment',x=px+.055,y=0,xend=px+.055,yend=py,linewidth=.3,colour=ink,
          arrow=grid::arrow(length=grid::unit(1.1,'mm'),ends='both'))+
 annotate('segment',x=px,y=py,xend=px+.07,yend=py,linewidth=.25,linetype=2,colour=grey)+
 annotate('text',x=px+.085,y=py/2,label='g',family='Arial',size=3.1,colour=ink)+
 annotate('segment',x=-.39,y=0,xend=-.39,yend=h,linewidth=.35,colour=ink,
          arrow=grid::arrow(length=grid::unit(1.3,'mm'),ends='both'))+
 annotate('text',x=-.42,y=h/2,label='h',family='Arial',size=3.2,colour=ink)+
 annotate('text',x=px/2+.025,y=(h+py)/2+.018,label='R',family='Arial',size=3.2,colour=ink)+
 annotate('text',x=.15,y=-.055,label='rigid plane',family='Arial',size=2.8,colour=ink)+
 coord_fixed(xlim=c(-.49,.43),ylim=c(-.07,.64),clip='off')+theme_void(base_family='Arial')+
 theme(plot.margin=margin(5,6,5,5))

# Figure 1b: assumed traction law, with no empirical interpretation.
g<-seq(-.006,.13,length.out=501)
traction<-ifelse(g<=0,-k*g,ifelse(g<a,-4*w*g/a^2*(1-(g/a)^2),0))
pb<-ggplot(data.frame(g,traction),aes(g,traction))+
 annotate('rect',xmin=0,xmax=a,ymin=-Inf,ymax=Inf,fill='#F8F2E4')+
 annotate('text',x=.05,y=.105,label='adhesive branch',family='Arial',size=3.0,colour=ink)+
 geom_hline(yintercept=0,linewidth=.24,colour=grey)+geom_vline(xintercept=c(0,a),linewidth=.3,linetype=2,colour=grey)+
 geom_line(linewidth=.75,colour=navy)+scale_x_continuous(breaks=c(0,.05,.1))+
 labs(x='Gap, g (µm)',y=expression('Traction, ' * -psi*"'" * ' (pN ' * mu * 'm'^{-2} * ')'))+
 theme_clean()+theme(legend.position='none')

# Figure 1c: model gap profiles for two prescribed proximity cases.
th<-seq(0,pi,length.out=501)
profiles<-rbind(data.frame(theta=th,g=R*cos(.15)+R*cos(th),case='Near g = 0'),
 data.frame(theta=th,g=a+R*cos(.15)+R*cos(th),case='Near g = a'))
pc<-ggplot(profiles,aes(theta,g,colour=case))+geom_hline(yintercept=c(0,a),linetype=2,colour=grey,linewidth=.3)+
 geom_line(linewidth=.7)+scale_colour_manual(values=c('Near g = 0'=navy,'Near g = a'=gold))+
 scale_x_continuous(breaks=c(0,pi/2,pi),labels=c('0','π/2','π'))+
 labs(x=expression('Angle, '*theta),y='Gap, g (µm)')+theme_clean()+theme(legend.position='bottom')
fig1<-plot_grid(pa,pb,pc,nrow=1,labels=c('a','b','c'),label_size=13,label_fontfamily='Arial',rel_widths=c(1,.99,1.04))

d<-read.csv(file.path(input,'45_indicador_gk_resultados.csv'),stringsAsFactors=FALSE)
cform<-read.csv(file.path(input,'48_forma_fechada_resultados.csv'),stringsAsFactors=FALSE)
mapping<-c(gk='Adaptive GK',guarded='Guarded GK',all='Split GK')
d$method<-factor(unname(mapping[d$method]),levels=names(pal))
base<-d[d$method=='Adaptive GK',]; floor<-1e-17

# Figure 2a: estimated error versus error relative to the reference.
base$actual_plot<-pmax(base$abs_error,floor);base$estimate_plot<-pmax(base$estimate,floor)
base$flag<-ifelse(tolower(as.character(base$false_safe))=='true','False safe',ifelse(tolower(as.character(base$warning))=='true','Warning','Other'))
qa<-ggplot(base,aes(estimate_plot,actual_plot))+
 annotate('rect',xmin=floor,xmax=tau,ymin=tau,ymax=1e-6,fill='#F8F2E4')+
 annotate('text',x=1e-16,y=1.2e-7,label='underestimated\nabove tolerance',family='Arial',size=3.0,hjust=0,lineheight=.9,colour=ink)+
 geom_abline(slope=1,intercept=0,colour=grey,linewidth=.28,linetype=2)+
 geom_vline(xintercept=tau,linewidth=.3,linetype=2,colour=ink)+geom_hline(yintercept=tau,linewidth=.3,linetype=2,colour=ink)+
 geom_point(aes(shape=flag,fill=flag),size=2.1,stroke=.35,colour=ink)+
 scale_shape_manual(values=c('False safe'=24,'Warning'=22,'Other'=21))+
 scale_fill_manual(values=c('False safe'=gold,'Warning'=blue,'Other'=navy))+
 scale_x_log10(limits=c(floor,1e-5),breaks=c(1e-16,1e-12,1e-8,1e-5))+
 scale_y_log10(limits=c(floor,1e-6),breaks=c(1e-16,1e-12,1e-8))+
 labs(x='GK estimated error (pN)',y='Error vs reference (pN)')+theme_clean()+
 theme(legend.position='bottom',legend.text=element_text(size=8.5))

# Figure 2b: all per-case errors, with deterministic offsets rather than stochastic jitter.
points<-data.frame(method=d$method,error=pmax(d$abs_error,floor),index=rep(seq_len(52),each=3))
points<-rbind(points,data.frame(method=factor('Closed form',levels=names(pal)),error=pmax(cform$abs_error,floor),index=seq_len(52)))
points$x<-as.numeric(points$method)+(points$index%%9-4)*.045
qb<-ggplot(points,aes(x,error,colour=method))+geom_hline(yintercept=tau,linetype=2,colour=ink,linewidth=.3)+
 annotate('text',x=2.55,y=3e-8,label='tolerance',family='Arial',size=3.0,hjust=0,colour=ink)+
 geom_point(size=1.25,alpha=.8)+scale_colour_manual(values=pal,drop=FALSE)+
 scale_x_continuous(breaks=1:4,labels=c('Adaptive\nGK','Guarded\nGK','Split\nGK','Closed\nform'))+
 scale_y_log10(limits=c(floor,1e-6),breaks=c(1e-16,1e-12,1e-8))+
 labs(x=NULL,y='Error vs reference (pN)')+theme_clean()+theme(legend.position='none')

# Figure 2c: explicit law-evaluation count (closed form has different operation types).
evals<-aggregate(evaluations~method,d,sum);evals$method<-factor(evals$method,levels=rev(names(pal)[1:3]))
qc<-ggplot(evals,aes(evaluations,method,fill=method))+geom_col(width=.58)+
 geom_text(aes(label=format(evaluations,big.mark=',')),hjust=-.12,size=3.15,family='Arial',colour=ink)+
 scale_fill_manual(values=pal,guide='none')+scale_x_continuous(limits=c(0,39000),expand=c(0,0))+
 labs(x='Contact-law evaluations',y=NULL)+theme_clean()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())

# Figure 2d: separate wall-only inclined-cylinder test from the earlier fixed protocol.
axial<-read.csv(file.path(input,'42_geometria_resultados.csv'),stringsAsFactors=FALSE)
axial<-subset(axial,abs(beta-.08)<1e-12)
axmax<-aggregate(force_error~method+n,axial,max)
axmax$method<-factor(axmax$method,levels=c('global','split'),labels=c('Global axial','Split axial'))
qd<-ggplot(axmax,aes(n,force_error,colour=method,shape=method))+
 geom_hline(yintercept=1.11347e-4,linetype=2,linewidth=.28,colour=ink)+
 geom_line(linewidth=.55)+geom_point(size=2)+
 scale_colour_manual(values=c('Global axial'=navy,'Split axial'=gold))+
 scale_shape_manual(values=c('Global axial'=16,'Split axial'=17))+
 scale_x_continuous(trans='log2',breaks=c(8,16,32,64))+
 scale_y_log10(breaks=c(1e-8,1e-6,1e-4))+
 labs(x='Axial Gauss points',y='Max error vs reference (pN)')+theme_clean()+
 theme(legend.position='bottom',legend.text=element_text(size=8.5))
fig2<-plot_grid(qa,qb,qc,qd,ncol=2,labels=c('a','b','c','d'),label_size=13,label_fontfamily='Arial',align='hv',axis='tblr')

# Figure 3a: constructed noncircular cross-section; translation h fixed by a sampled case.
nc<-read.csv(file.path(input,'51_nao_circular_resultados.csv'),stringsAsFactors=FALSE)
sample_case<-subset(nc,id=='b0.0_pi_e0.13' & method=='gk')
stopifnot(nrow(sample_case)==1)
theta_nc<-seq(0,2*pi,length.out=501); rx<-.36;ry<-.28209479177387814;q<-.045
contour<-data.frame(x=rx*sin(theta_nc),y=sample_case$h+ry*cos(theta_nc)+q*cos(2*theta_nc))
ra<-ggplot()+annotate('rect',xmin=-.43,xmax=.43,ymin=-.035,ymax=0,fill='#EFF1F4')+
 geom_hline(yintercept=0,linewidth=.5,colour=ink)+geom_path(data=contour,aes(x,y),colour=navy,linewidth=.7)+
 coord_fixed(xlim=c(-.43,.43),ylim=c(-.07,.61),clip='off')+
 annotate('text',x=.19,y=-.055,label='rigid plane',family='Arial',size=2.8,colour=ink)+
 theme_void(base_family='Arial')+theme(plot.margin=margin(5,6,5,5))

# Figure 3b: paired per-case error across numerical methods.
nc$method_label<-factor(unname(mapping[nc$method]),levels=c('Adaptive GK','Guarded GK','Split GK'))
nc$error_plot<-pmax(nc$abs_error,floor)
rb<-ggplot(nc,aes(method_label,error_plot,group=id))+
 geom_hline(yintercept=tau,linetype=2,colour=ink,linewidth=.3)+
 annotate('text',x='Adaptive GK',y=2e-8,label='tolerance',family='Arial',size=3.0,hjust=0,colour=ink)+
 geom_line(colour=grey,alpha=.28,linewidth=.22)+
 geom_point(aes(colour=method_label),size=1.55,alpha=.9)+
 scale_colour_manual(values=pal,guide='none')+
 scale_y_log10(limits=c(floor,1e-7),breaks=c(1e-16,1e-12,1e-8))+
 labs(x=NULL,y='Error vs reference (pN)')+theme_clean()+
 theme(axis.text.x=element_text(size=9))

# Figure 3c: both integrand geometry and additional root-search geometry calls.
ncost<-aggregate(cbind(law_evals,root_geometry_evals)~method_label,nc,sum)
ncost$method_label<-factor(ncost$method_label,levels=rev(c('Adaptive GK','Guarded GK','Split GK')))
costlong<-rbind(data.frame(method=ncost$method_label,part='During force integration',calls=ncost$law_evals),
 data.frame(method=ncost$method_label,part='Root search',calls=ncost$root_geometry_evals))
costlong$part<-factor(costlong$part,levels=c('Root search','During force integration'))
ncost$total<-ncost$law_evals+ncost$root_geometry_evals
rc<-ggplot(costlong,aes(calls,method,fill=part))+geom_col(width=.58,position='stack')+
 geom_text(data=ncost,aes(x=total,y=method_label,label=format(total,big.mark=',')),inherit.aes=FALSE,
           hjust=-.15,size=3.15,family='Arial',colour=ink)+
 scale_fill_manual(values=c('During force integration'=navy,'Root search'=gold))+
 scale_x_continuous(limits=c(0,14000),expand=c(0,0))+
 labs(x='Geometry evaluations across 18 cases',y=NULL)+theme_clean()+
 theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),legend.position='bottom')
fig3top<-plot_grid(ra,rb,nrow=1,labels=c('a','b'),label_size=13,label_fontfamily='Arial',rel_widths=c(.85,1.15))
fig3bottom<-plot_grid(rc,nrow=1,labels='c',label_size=13,label_fontfamily='Arial')
fig3<-plot_grid(fig3top,fig3bottom,ncol=1,rel_heights=c(.58,.42))

# Figure 4: post-review diagnostics are exploratory, not the primary protocol.
review_dir <- file.path(paper, 'review_checks')
diag <- read.csv(file.path(review_dir,'isolated_breakpoints.csv'),stringsAsFactors=FALSE)
case_ids <- unique(diag$id)
failure_rows <- d[as.character(d$method)=='Adaptive GK' & d$id %in% case_ids,]
stopifnot(length(case_ids)==2, nrow(failure_rows)==2)
first <- failure_rows[grepl('_e0.03$',failure_rows$id),]
theta_plot <- seq(0,pi,length.out=601)
gap_plot <- data.frame(theta=theta_plot,g=first$h+first$r*cos(theta_plot))
t_a <- acos((a-first$h)/first$r); t_0 <- acos(-first$h/first$r)
sa <- ggplot(gap_plot,aes(theta,g))+
 annotate('rect',xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=0,fill='#F8F2E4')+
 geom_hline(yintercept=c(0,a),linetype=2,linewidth=.3,colour=grey)+
 geom_line(linewidth=.8,colour=navy)+
 annotate('point',x=t_a,y=a,shape=24,fill=gold,colour=ink,size=2.5)+
 annotate('point',x=t_0,y=0,shape=21,fill=blue,colour=ink,size=2.5)+
 annotate('text',x=.85,y=-.40,label='g < 0: penalty overlap',size=2.8,family='Arial',colour=ink)+
 annotate('text',x=.25,y=.16,label='g = a',size=2.8,family='Arial',colour=ink)+
 annotate('text',x=1.2,y=.055,label='g = 0',size=2.8,family='Arial',colour=ink)+
 scale_x_continuous(breaks=c(0,pi/2,pi),labels=c('0','π/2','π'))+
 labs(x=expression('Angle, '*theta),y='Gap, g (µm)')+theme_clean()+theme(legend.position='none')
orig <- data.frame(id=failure_rows$id,strategy='none',abs_error_pN=failure_rows$abs_error)
iso <- rbind(orig,diag[,c('id','strategy','abs_error_pN')])
iso$strategy <- factor(iso$strategy,levels=c('none','near_adhesive_only','interior_compression_only','both'),labels=c('None','g = a','g = 0','Both'))
iso$case <- ifelse(grepl('_e0.03$',iso$id),'0.03 rad','0.15 rad')
iso$error_plot <- pmax(iso$abs_error_pN,1e-17)
sb <- ggplot(iso,aes(strategy,error_plot,group=case,colour=case,shape=case))+
 geom_hline(yintercept=tau,linetype=2,colour=ink,linewidth=.3)+
 geom_line(linewidth=.5,alpha=.8)+geom_point(size=2.2)+
 scale_colour_manual(values=c('0.03 rad'=navy,'0.15 rad'=gold))+
 scale_shape_manual(values=c('0.03 rad'=16,'0.15 rad'=17))+
 scale_y_log10(limits=c(1e-17,1e-6),breaks=c(1e-16,1e-12,1e-8))+
 labs(x='Supplied breakpoints',y='Error vs reference (pN)')+theme_clean()+
 theme(legend.text=element_text(size=8.5))
fixed <- read.csv(file.path(review_dir,'fixed_segment_gauss.csv'),stringsAsFactors=FALSE)
make_fixed_panel <- function(which_geometry, primary) {
 df <- subset(fixed,geometry==which_geometry)
 # Sum and maximum are calculated directly to keep the units of each metric explicit.
 rows <- do.call(rbind,lapply(c(4,8),function(order) {
   x <- subset(df,order_per_segment==order);
   data.frame(method=paste('Fixed',order),error=max(x$abs_error_pN),
              count=sum(x$force_evaluations),fail=sum(x$abs_error_pN>tau),n=nrow(x))
 }))
 rows <- rbind(rows,data.frame(method='Split GK',error=max(primary$abs_error),
   count=sum(if ('evaluations' %in% names(primary)) primary$evaluations else primary$law_evals),
   fail=sum(primary$abs_error>tau),n=nrow(primary)))
 rows$method <- factor(rows$method,levels=c('Fixed 4','Fixed 8','Split GK'))
 rows$label <- paste0(rows$fail,'/',rows$n)
 rows$error <- pmax(rows$error,1e-17)
 g <- ggplot(rows,aes(method,error))+geom_hline(yintercept=tau,linetype=2,colour=ink,linewidth=.3)+
   geom_point(aes(colour=method,shape=method),size=2.7)+
   geom_text(aes(label=label),vjust=-.8,size=2.9,family='Arial',colour=ink)+
   scale_colour_manual(values=c('Fixed 4'=grey,'Fixed 8'=gold,'Split GK'=navy),guide='none')+
   scale_shape_manual(values=c('Fixed 4'=15,'Fixed 8'=17,'Split GK'=16),guide='none')+
   scale_y_log10(limits=if(which_geometry=='circular') c(1e-17,1e-3) else c(1e-17,1e1),
                 breaks=if(which_geometry=='circular') c(1e-16,1e-12,1e-8,1e-4) else c(1e-16,1e-12,1e-8,1e-4,1))+
   labs(x='Integration rule',y='Max error vs reference (pN)')+theme_clean()+
   theme(axis.text.x=element_text(size=8.5),legend.position='none',plot.margin=margin(12,5,5,5))
 g
}
sc <- make_fixed_panel('circular',d[as.character(d$method)=='Split GK',])
sd <- make_fixed_panel('noncircular',nc[as.character(nc$method_label)=='Split GK',])
fig4<-plot_grid(sa,sb,sc,sd,ncol=2,labels=c('a','b','c','d'),label_size=13,label_fontfamily='Arial',align='hv',axis='tblr')

save_figure<-function(plot,name,width,height){
 ggsave(file.path(here,paste0(name,'.png')),plot=plot,width=width,height=height,units='mm',dpi=450,bg='white',limitsize=FALSE)
 ggsave(file.path(here,paste0(name,'.svg')),plot=plot,width=width,height=height,units='mm',device=svglite::svglite,bg='white',limitsize=FALSE)
 ggsave(file.path(here,paste0(name,'.pdf')),plot=plot,width=width,height=height,units='mm',device=cairo_pdf,bg='white',limitsize=FALSE)
}
save_figure(fig1,'figure1',180,65)
save_figure(fig2,'figure2',180,150)
save_figure(fig3,'figure3',180,140)
save_figure(fig4,'figure4',180,150)

# Separately labelled vector panels provide an editable-artwork option.
panel_dir <- file.path(here, 'submission_panels')
dir.create(panel_dir, showWarnings=FALSE)
panel_export <- function(plot, name, label, width, height) {
  labelled <- plot_grid(plot, labels=label, label_size=13, label_fontfamily='Arial')
  ggsave(file.path(panel_dir,paste0(name,'.pdf')), plot=labelled, width=width, height=height, units='mm', device=cairo_pdf, bg='white')
  ggsave(file.path(panel_dir,paste0(name,'.svg')), plot=labelled, width=width, height=height, units='mm', device=svglite::svglite, bg='white')
}
panel_export(pa,'Figure_1a','a',65,65);panel_export(pb,'Figure_1b','b',65,65);panel_export(pc,'Figure_1c','c',65,65)
panel_export(qa,'Figure_2a','a',90,75);panel_export(qb,'Figure_2b','b',90,75);panel_export(qc,'Figure_2c','c',90,75);panel_export(qd,'Figure_2d','d',90,75)
panel_export(ra,'Figure_3a','a',80,80);panel_export(rb,'Figure_3b','b',100,80);panel_export(rc,'Figure_3c','c',180,60)

panel_export(sa,'Figure_4a','a',90,75);panel_export(sb,'Figure_4b','b',90,75);panel_export(sc,'Figure_4c','c',90,75);panel_export(sd,'Figure_4d','d',90,75)
