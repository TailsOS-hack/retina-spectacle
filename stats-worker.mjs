import {compareGenes} from './research.mjs';
self.onmessage=({data})=>{try {self.postMessage({id:data.id,rows:compareGenes(data.geneNames,data.vectors,data.a,data.b)});}catch(error){self.postMessage({id:data.id,error:error.message});}};
