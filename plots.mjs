import {palette,esc,fmt,summarize} from './research.mjs';
const FONT='12px -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif';
export function colorValue(value,max) {if(value<=0)return '#e4e9e6';const t=Math.max(0,Math.min(1,value/(max||1)));return `rgb(${Math.round(229-205*t)},${Math.round(240-128*t)},${Math.round(233-137*t)})`;}
export function setupCanvas(canvas,width,height) {
  const dpr=Math.min(globalThis.devicePixelRatio||1,2);canvas.style.width=width+'px';canvas.style.height=height+'px';canvas.width=Math.round(width*dpr);canvas.height=Math.round(height*dpr);
  const ctx=canvas.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);ctx.font=FONT;ctx.fillStyle='#fff';ctx.fillRect(0,0,width,height);return ctx;
}
function wrap(ctx,text,x,y,width,lineHeight=16,maxLines=3) {
  const words=String(text).split(/\s+/);let line='',row=0;
  for(const word of words){const next=(line+' '+word).trim();if(ctx.measureText(next).width>width&&line){ctx.fillText(line,x,y+row*lineHeight);line=word;row++;if(row>=maxLines-1)break;}else line=next;}
  ctx.fillText(line,x,y+row*lineHeight);
}
export function drawDistribution(canvas,groups,values,{width=850,kind='box',showPoints=true,gene='',scale='ln(1 + CP10K)'}={}) {
  const widthActual=Math.max(width,groups.length*140+110),h=440,L=68,R=widthActual-28,T=48,B=332,ctx=setupCanvas(canvas,widthActual,h),step=(R-L)/Math.max(1,groups.length);
  const summaries=groups.map(g=>summarize(g.ids,values)),maximum=Math.max(1,...summaries.map(s=>s.max||0))*1.08,y=v=>B-v/maximum*(B-T),points=[];
  ctx.textAlign='left';ctx.fillStyle='#263f36';ctx.font='600 14px -apple-system, sans-serif';ctx.fillText(gene,L,24);ctx.font=FONT;ctx.fillStyle='#69776f';ctx.fillText(scale,Math.max(L+110,R-180),24);
  for(let k=0;k<=5;k++){const v=maximum*k/5,yy=y(v);ctx.strokeStyle='#e9eeeb';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(L,yy);ctx.lineTo(R,yy);ctx.stroke();ctx.fillStyle='#69776f';ctx.textAlign='right';ctx.fillText(v.toFixed(1),L-12,yy+4);}
  groups.forEach((g,j)=>{
    const s=summaries[j],x=L+step*(j+.5),color=palette[j%palette.length];ctx.textAlign='center';ctx.fillStyle='#34483f';wrap(ctx,g.name,x,B+26,step-14);ctx.fillStyle='#7b857f';ctx.fillText(`n = ${fmt(g.ids.length)}`,x,B+82);
    if(!s.n){ctx.fillText('No cells',x,(T+B)/2);return;}
    if(kind==='violin'){
      const bins=new Float64Array(40);for(const i of g.ids)bins[Math.min(39,Math.floor(values[i]/maximum*40))]++;
      const smooth=bins.map((_,i)=>{let v=0;for(let k=-3;k<=3;k++)if(i+k>=0&&i+k<40)v+=bins[i+k]*Math.exp(-k*k/3);return v;}),max=Math.max(...smooth,1),half=Math.min(44,step*.25);
      ctx.beginPath();for(let i=0;i<40;i++){const xx=x-smooth[i]/max*half,yy=y((i+.5)/40*maximum);i?ctx.lineTo(xx,yy):ctx.moveTo(xx,yy);}for(let i=39;i>=0;i--)ctx.lineTo(x+smooth[i]/max*half,y((i+.5)/40*maximum));ctx.closePath();ctx.fillStyle=color+'25';ctx.fill();ctx.strokeStyle=color;ctx.stroke();
    }
    if(showPoints){ctx.fillStyle=color;ctx.globalAlpha=Math.min(.6,Math.max(.28,500/Math.max(500,g.ids.length)));for(const i of g.ids){const hash=((Math.imul(i+1,2654435761)>>>0)%10001)/10000,xx=x+(hash-.5)*Math.min(step*.68,104),yy=y(values[i]);ctx.fillRect(xx-1,yy-1,2,2);points.push([xx,yy,i,j]);}ctx.globalAlpha=1;}
    const hw=kind==='box'?Math.min(30,step*.19):6;ctx.strokeStyle=color;ctx.lineWidth=1.6;ctx.beginPath();ctx.moveTo(x,y(s.low));ctx.lineTo(x,y(s.high));ctx.moveTo(x-hw*.5,y(s.low));ctx.lineTo(x+hw*.5,y(s.low));ctx.moveTo(x-hw*.5,y(s.high));ctx.lineTo(x+hw*.5,y(s.high));ctx.stroke();ctx.fillStyle=color+'24';ctx.fillRect(x-hw,y(s.q3),hw*2,Math.max(1,y(s.q1)-y(s.q3)));ctx.strokeRect(x-hw,y(s.q3),hw*2,Math.max(1,y(s.q1)-y(s.q3)));ctx.lineWidth=2.6;ctx.beginPath();ctx.moveTo(x-hw,y(s.median));ctx.lineTo(x+hw,y(s.median));ctx.stroke();
  });
  if(!groups.length){ctx.textAlign='center';ctx.fillStyle='#718078';ctx.fillText('No groups selected. Include a group or reset the filters.',widthActual/2,h/2);}
  canvas.dataset.cellCount=String(groups.reduce((n,g)=>n+g.ids.length,0));canvas.dataset.pointsDrawn=String(points.length);canvas.dataset.groupCount=String(groups.length);
  return {points,summaries};
}
export function panelSVG(groups,vectors,geneNames,kind='dot') {
  const L=245,step=92,row=45,w=Math.max(670,L+step*geneNames.length+24),h=115+groups.length*row,grid=groups.map(g=>vectors.map(v=>summarize(g.ids,v))),max=Math.max(.01,...grid.flat().map(s=>s.mean||0));
  let svg=`<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" role="img" aria-label="${kind==='dot'?'Dot plot':'Heatmap'} of selected genes and groups" style="font-family:-apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;background:white"><rect width="100%" height="100%" fill="white"/>`;
  geneNames.forEach((name,j)=>svg+=`<text x="${L+step*(j+.5)}" y="28" text-anchor="middle" font-size="12" font-weight="600" fill="#294137">${esc(name)}</text>`);
  groups.forEach((g,i)=>{
    const cy=58+i*row;svg+=`<rect x="0" y="${cy-19}" width="${w}" height="${row}" fill="${i%2?'#f8faf8':'#fff'}"/><text x="12" y="${cy+4}" font-size="12" fill="#43554b"><title>${esc(g.name)} · ${fmt(g.ids.length)} cells</title>${esc(g.name.length>31?g.name.slice(0,29)+'…':g.name)}</text>`;
    grid[i].forEach((s,j)=>{const title=`${g.name} / ${geneNames[j]}: n=${s.n}, mean=${s.mean?.toFixed(4)??'unavailable'}, detected=${s.pct?.toFixed(2)??'unavailable'}%`,x=L+step*(j+.5);if(kind==='heatmap')svg+=`<rect x="${x-step/2+3}" y="${cy-16}" width="${step-6}" height="32" fill="${s.n?colorValue(s.mean,max):'#f1f1f1'}"><title>${esc(title)}</title></rect>`;else svg+=`<circle cx="${x}" cy="${cy}" r="${s.n?Math.sqrt(s.pct/100)*15:0}" fill="${colorValue(s.mean,max)}"><title>${esc(title)}</title></circle>${!s.n?`<text x="${x}" y="${cy+4}" text-anchor="middle" fill="#aaa">—</text>`:''}`;});
  });
  const y=h-26;svg+=`<text x="12" y="${y}" font-size="11" fill="#6b7b71">Mean expression: 0</text>`;for(let i=0;i<70;i++)svg+=`<rect x="${125+i}" y="${y-10}" width="1.2" height="10" fill="${colorValue(i/69*max,max)}"/>`;svg+=`<text x="203" y="${y}" font-size="11" fill="#6b7b71">${max.toFixed(2)}</text>`;
  if(kind==='dot'){svg+=`<text x="285" y="${y}" font-size="11" fill="#6b7b71">Detected cells</text>`;[25,50,100].forEach((pct,i)=>svg+=`<circle cx="${410+i*83}" cy="${y-6}" r="${Math.sqrt(pct/100)*12}" fill="#568c78"/><text x="${427+i*83}" y="${y}" font-size="11" fill="#6b7b71">${pct}%</text>`);}
  return {svg:svg+'</svg>',grid};
}
export function drawMap(canvas,meta,cells,ids,values,{width=950,height=530,projection='umap',color='condition',sampleInfo=[],typeNames=new Map(),selectionA=new Set(),selectionB=new Set(),pan={x:0,y:0,z:1},subset=null,path=[]}={}) {
  const ctx=setupCanvas(canvas,width,height),L=45,R=width-30,T=25,B=height-43,ts=projection==='tsne',xv=subset?subset[ts?'tx':'x']:cells[ts?'tx':'x'],yv=subset?subset[ts?'ty':'y']:cells[ts?'ty':'y'],at=i=>subset?subset.lookup.get(i):i;
  let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity;for(let j=0;j<xv.length;j++){xmin=Math.min(xmin,xv[j]);xmax=Math.max(xmax,xv[j]);ymin=Math.min(ymin,yv[j]);ymax=Math.max(ymax,yv[j]);}
  const scale=Math.min((R-L-20)/Math.max(.01,xmax-xmin),(B-T-20)/Math.max(.01,ymax-ymin))*pan.z,mx=(L+R)/2+pan.x,my=(T+B)/2+pan.y;
  ctx.strokeStyle='#e8ede9';ctx.lineWidth=1;for(let k=0;k<5;k++){const x=L+(R-L)*k/4,y=T+(B-T)*k/4;ctx.beginPath();ctx.moveTo(x,T);ctx.lineTo(x,B);ctx.moveTo(L,y);ctx.lineTo(R,y);ctx.stroke();}
  const labels=[];const getLabel=i=>color==='subset'&&subset?`Reclustered C${subset.cluster[at(i)]}`:color==='sample'?sampleInfo[cells.sample[i]].title:color==='celltype'?typeNames.get(cells.cluster[i]):color==='cluster'?`C${cells.cluster[i]}`:sampleInfo[cells.sample[i]].condition;
  for(const i of ids){const label=getLabel(i);if(!labels.includes(label))labels.push(label);}
  let max=0;for(const i of ids)max=Math.max(max,values?.[i]||0);
  const order=color==='gene'&&values?[...ids].sort((a,b)=>values[a]-values[b]):ids,points=[];ctx.save();ctx.beginPath();ctx.rect(L,T,R-L,B-T);ctx.clip();ctx.globalAlpha=.76;
  for(const i of order){const j=at(i);if(j===undefined)continue;const x=(xv[j]-(xmin+xmax)/2)*scale+mx,y=my-(yv[j]-(ymin+ymax)/2)*scale;if(x<L||x>R||y<T||y>B)continue;ctx.fillStyle=color==='gene'?colorValue(values?.[i]||0,max):color==='selection'?(selectionA.has(i)?palette[0]:selectionB.has(i)?palette[1]:'#dbe1dd'):palette[labels.indexOf(getLabel(i))%palette.length];const r=Math.max(1.2,Math.min(2.8,1.5*Math.sqrt(pan.z)));ctx.fillRect(x-r,y-r,r*2,r*2);points.push([x,y,i]);}
  ctx.globalAlpha=1;if(path.length){ctx.strokeStyle='#276951';ctx.lineWidth=2;ctx.setLineDash([4,3]);ctx.beginPath();path.forEach(([x,y],i)=>i?ctx.lineTo(x,y):ctx.moveTo(x,y));ctx.stroke();ctx.setLineDash([]);}ctx.restore();
  ctx.fillStyle='#748276';ctx.textAlign='center';ctx.fillText(projection.toUpperCase()+' 1',(L+R)/2,height-12);ctx.save();ctx.translate(13,(T+B)/2);ctx.rotate(-Math.PI/2);ctx.fillText(projection.toUpperCase()+' 2',0,0);ctx.restore();
  if(!ids.length){ctx.fillStyle='#69776f';ctx.fillText('No cells match your filters.',width/2,height/2);}
  return {points,labels,max};
}
