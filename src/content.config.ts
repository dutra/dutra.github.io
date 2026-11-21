import { defineCollection } from 'astro:content';
import { glob } from 'astro/loaders';
import { z } from 'astro/zod';

const publications = defineCollection({
  loader: glob({ pattern: '**/*.md', base: './src/content/publications' }),
  schema: z.object({
    title: z.string(),
    /** Comma-separated author list, in order. */
    author: z.string(),
    /** Year of publication (or of the preprint). */
    date: z.string(),
    journal: z.string(),
    arxiv: z.string().optional(),
    doi: z.string().optional(),
    first_author: z.boolean().default(false),
    /** Short summary shown under "Abstract". */
    description: z.string(),
  }),
});

export const collections = { publications };
